package es.elixir.bsc.ngs.nova.gecoz.tools;

import es.elixir.bsc.ngs.nova.algo.ssa.GSSA;
import es.elixir.bsc.ngs.nova.gecoz.GecozFileReader;
import es.elixir.bsc.ngs.nova.gecoz.GecozRefBlockHeader;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import static java.nio.file.StandardOpenOption.CREATE;
import static java.nio.file.StandardOpenOption.READ;
import static java.nio.file.StandardOpenOption.TRUNCATE_EXISTING;
import static java.nio.file.StandardOpenOption.WRITE;
import java.util.EnumSet;
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
            
            FileChannel channel = FileChannel.open(opath, EnumSet.of(CREATE,READ,WRITE, TRUNCATE_EXISTING));
            ByteBuffer buf = channel.map(FileChannel.MapMode.READ_WRITE, 0, Math.min(to, bheader.len - 1) - from);

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
    
    static void fasta(Path ipath, Path opath) {
        
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

            try(BufferedOutputStream out = new BufferedOutputStream(Files.newOutputStream(opath))) {
                final long t1 = System.nanoTime();

                for (GecozRefBlockHeader bheader : reader.getBlockHeaders()) {
                    GSSA ssa = reader.read(bheader);
                    if (ssa == null) {
                        System.err.println("no block found " + bheader.headers[0] + " skipping ...");
                        continue;
                    }

                    for (String header : bheader.headers) {

                        System.out.print("writing '>" + header + "'...");
                        out.write('>');
                        out.write(header.getBytes("UTF8"));
                        out.write('\n');

                        final int nstr = bheader.findHeader(header);
                        
                        if (nstr < 0) {
                            System.err.println("no sequence found: " + header + " skipping...");
                            continue;
                        }
                        
                        long from = 0;

                        byte[] arr = new byte[1024 * 1024 * 4]; // 4Mb buffer
                        ByteBuffer buf = ByteBuffer.wrap(arr);
                        do {
                            buf.rewind();
                            ssa.extract(buf, nstr, from);
                            for (int i = 0, n = buf.position() & 0xFFFFFFC0; i < n; i += 64) {
                                out.write(arr, i, 64);
                                out.write('\n');
                            }
                            from += buf.position();
                        } while (buf.position() == buf.limit());
                        final int tail = buf.position() & 63;
                        if (tail > 0) {
                            out.write(arr, buf.position() - tail, tail);
                            out.write('\n');                                
                        }
                        System.out.println(" ok.");
                    }
                }
                final long t2 = System.nanoTime();
                System.out.println("finished in " + ((t2 - t1)/1000000) + " ms.");
            } catch (IOException ex) {
                System.err.println("error extracting fasta to " + opath);
                System.exit(1);                        
            }
        } catch(IOException | DataFormatException ex) {
            System.err.println("error reading a file: " + ipath);
            System.exit(1);
        }
    }
}
