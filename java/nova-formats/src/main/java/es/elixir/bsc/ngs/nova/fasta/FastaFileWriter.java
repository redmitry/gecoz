package es.elixir.bsc.ngs.nova.fasta;

import java.io.Closeable;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.Pipe;
import java.nio.channels.Pipe.SinkChannel;
import java.nio.channels.Pipe.SourceChannel;
import java.nio.file.Path;
import static java.nio.file.StandardOpenOption.CREATE;
import static java.nio.file.StandardOpenOption.READ;
import static java.nio.file.StandardOpenOption.TRUNCATE_EXISTING;
import static java.nio.file.StandardOpenOption.WRITE;
import java.util.EnumSet;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Dmitry Repchevsky
 */

public class FastaFileWriter implements Closeable {
    
    public final static int LINE_LENGTH = 50;
    
    private final ExecutorService executor;
    private final FileChannel channel;

    public FastaFileWriter(Path path) throws IOException {
        this(path, 1);
    }
    
    public FastaFileWriter(Path path, int threads) throws IOException {
        channel = FileChannel.open(path, EnumSet.of(CREATE,READ,WRITE,TRUNCATE_EXISTING));
        
        threads = Math.min(threads, Runtime.getRuntime().availableProcessors());
        
        Logger.getLogger(FastaFileWriter.class.getName()).log(Level.FINER, "fasta writer use {0} threads\n", threads);

        executor = new ThreadPoolExecutor(1, threads, 0L, TimeUnit.MILLISECONDS, 
                                          new ArrayBlockingQueue<>(1), 
                                          new RejectedExecutionHandler() {
            @Override
            public void rejectedExecution(Runnable runnable, ThreadPoolExecutor executor) {
                if (!executor.isShutdown()) {
                    try {
                        executor.getQueue().put(runnable);
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                        Logger.getLogger(FastaFileWriter.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
        });
    }
    
    public SinkChannel write(TFastaSequence seq) throws IOException {
        
        Logger.getLogger(FastaFileWriter.class.getName()).log(Level.FINER, "writing {0} ({1} bytes)\n", new Object[]{seq.header, Integer.toString(seq.length)});
        
        final byte[] header = seq.header.getBytes("UTF8");
        final ByteBuffer buf = ByteBuffer.allocate(header.length + 2).put((byte)'>').put(header).put((byte)'\n');
        buf.rewind();
        while(channel.write(buf) > 0);

        final long pos = channel.position();
        final long size = seq.multiline ? seq.length + (seq.length / LINE_LENGTH) + 1 : seq.length + 1;

        final ByteBuffer out = channel.map(FileChannel.MapMode.READ_WRITE, pos, size);
        channel.position(pos + size);
                
        final Pipe pipe = Pipe.open();
        
        executor.submit(new FastaSequenceWriter(pipe.source(), out));

        return pipe.sink();
    }
    
    @Override
    public void close() throws IOException {
        try {
            executor.shutdown();
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch (InterruptedException ex) {
            throw new IOException(ex.getMessage());
        } finally {
            channel.close();
        }
    }

    public static class FastaSequenceWriter implements Runnable {

        private final SourceChannel src;
        private final ByteBuffer out;
        
        public FastaSequenceWriter(SourceChannel src, ByteBuffer out) {
            this.src = src;
            this.out = out;
            
            out.limit(LINE_LENGTH);
        }
        
        @Override
        public void run() {
            try {
                do {
                    while (out.remaining() > 1) { // last '\n'
                        src.read(out);
                    }
                    out.limit(Math.min(out.capacity(), out.limit() + LINE_LENGTH + 1));
                    out.put((byte)'\n');
                } while (out.remaining() > 0);
            } catch (Throwable ex) {
                Logger.getLogger(FastaFileWriter.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
}
