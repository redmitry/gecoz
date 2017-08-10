package es.elixir.bsc.ngs.nova.gecoz.tools;

import es.elixir.bsc.ngs.nova.algo.ssa.GSSA;
import es.elixir.bsc.ngs.nova.fasta.FastaFileWriter;
import es.elixir.bsc.ngs.nova.fasta.TFastaSequence;
import es.elixir.bsc.ngs.nova.gecoz.GecozFileReader;
import es.elixir.bsc.ngs.nova.gecoz.GecozRefBlockHeader;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.Pipe.SinkChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import static java.nio.file.StandardOpenOption.CREATE;
import static java.nio.file.StandardOpenOption.READ;
import static java.nio.file.StandardOpenOption.TRUNCATE_EXISTING;
import static java.nio.file.StandardOpenOption.WRITE;
import java.util.EnumSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;

/**
 * Auxiliary class to extract sequences from GecoZ files.
 * 
 * @author Dmitry Repchevsky
 */
public class GecoRead {

    static void sequence(Path ipath, String header, long from, long to, Path opath) {
        try {
            if (!Files.exists(ipath) || Files.isDirectory(ipath)) {
                System.err.println("no gecoz file found: " + ipath);
                System.exit(1); 
            }
            
            if (!GecozFileReader.checkFormat(ipath)) {
                System.err.println("invalid gecoz file format: " + ipath);
                System.exit(1);
            }

            GecozFileReader reader = new GecozFileReader(ipath);

            GecozRefBlockHeader bheader = reader.findBlockHeader(header);
            if (bheader == null) {
                System.err.println("no sequence found: " + header);
                System.exit(1);
            }

            GSSA ssa = reader.read(bheader);
            if (ssa == null) {
                System.err.println("no block found: " + bheader.len + " bytes");
                System.exit(1);
            }

            final int nstr = bheader.findHeader(header);
            if (nstr < 0) {
                System.err.println("no sequence found: " + header);
                System.exit(1);                
            }
            
            to = Math.min(to, ssa.getLength(nstr));

            FileChannel channel = FileChannel.open(opath, EnumSet.of(CREATE,READ,WRITE, TRUNCATE_EXISTING));
            ByteBuffer buf = channel.map(FileChannel.MapMode.READ_WRITE, 0, to - from);

            System.out.println("extracting '" + header + "' (from " + from + " to " + (to == Integer.MAX_VALUE ? ".." : Long.toString(to)) + ")");
            final long t1 = System.nanoTime();
            ssa.extract(buf, nstr, from);
            final long t2 = System.nanoTime();

            System.out.println("finished in " + ((t2 - t1)/1000000) + " ms.");

        } catch(IOException | DataFormatException ex) {
            System.err.println("error reading a file: " + ipath);
            System.exit(1);
        }
    }
    
    static void fasta(Path ipath, Path opath, int threads) {
        
        try {
            if (!Files.exists(ipath) || Files.isDirectory(ipath)) {
                System.err.println("no gecoz file found: " + ipath);
                System.exit(1); 
            }
            
            if (!GecozFileReader.checkFormat(ipath)) {
                System.err.println("invalid gecoz file format: " + ipath);
                System.exit(1);
            }

            ExecutorService executor = Executors.newFixedThreadPool(threads);
            
            final long t1 = System.nanoTime();
            try(GecozFileReader reader = new GecozFileReader(ipath);
                FastaFileWriter writer = new FastaFileWriter(opath, threads)) {
                for (GecozRefBlockHeader bheader : reader.getBlockHeaders()) {
                    GSSA ssa = reader.read(bheader);
                    if (ssa == null) {
                        System.err.println("no block found " + bheader.headers[0] + " skipping ...");
                        continue;
                    }

                    for (String header : bheader.headers) {
                        
                        final int nstr = bheader.findHeader(header);
                        if (nstr < 0) {
                            System.err.println("no sequence found: " + header + " skipping...");
                            continue;
                        }
                        
                        final int len = (int)ssa.getLength(nstr);
                        final SinkChannel sink = writer.write(new TFastaSequence(header, (int)len, true));
                        
                        executor.submit(new SequenceExtractor(ssa, nstr, len, sink));
                    }
                }
            } catch (IOException ex) {
                System.err.println("error extracting fasta to " + opath);
                System.exit(1);                        
            } finally {
                executor.shutdown();
                try {
                    executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
                } catch (InterruptedException ex) {
                    Logger.getLogger(GecoRead.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            final long t2 = System.nanoTime();
            System.out.println("finished in " + ((t2 - t1)/1000000) + " ms.");
        } catch(IOException | DataFormatException ex) {
            System.err.println("error reading a file: " + ipath);
            System.exit(1);
        }
    }
    
    public static class SequenceExtractor implements Runnable {

        private final GSSA ssa;
        private final int nstr;
        private final int len;
        private final SinkChannel sink;
        
        public SequenceExtractor(GSSA ssa, int nstr, int len, SinkChannel sink) {
            this.ssa = ssa;
            this.nstr = nstr;
            this.len = len;
            this.sink = sink;
        }
        
        @Override
        public void run() {
            
            try {
                long from = 0;

                ByteBuffer buf = ByteBuffer.allocate(1024 * 1024 * 4); // 4Mb buffer
                do {
                    buf.rewind();
                    ssa.extract(buf, nstr, from);
                    buf.limit(buf.position());
                    buf.rewind();
                    sink.write(buf);
                    from += buf.position();
                } while (from < len);
            }
            catch (IOException ex) {
                Logger.getLogger(GecoIndex.class.getName()).log(Level.SEVERE, ex.getMessage());
            }
        }
    }
}
