package es.elixir.bsc.ngs.nova.fasta;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
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
    
    private ExecutorService executor;
    private final FileChannel channel;

    public FastaFileWriter(Path path) throws IOException {
        this(path, 1);
    }
    
    public FastaFileWriter(Path path, int threads) throws IOException {
        channel = FileChannel.open(path, EnumSet.of(CREATE,READ,WRITE,TRUNCATE_EXISTING));
        
        threads = Math.min(threads, Runtime.getRuntime().availableProcessors());
        
        Logger.getLogger(FastaFileWriter.class.getName()).log(Level.FINER, "fasta writer use {0} threads\n", threads);

        if (threads > 1) {
            setExecutor(threads);
        }
    }

    /**
     * Directly writes the fasta sequence to the file.
     * 
     * @param header fasta header
     * @param seq sequence array
     * @param off sequence array offset
     * @param len length of sequence
     * 
     * @throws IOException 
     */
    public void write(final String header, final byte[] seq, final int off, final int len) throws IOException {
        write(header, false, seq, off, len);
    }
    
    private void write(final String header, final boolean multiline, final byte[] seq, final int off, final int len) throws IOException {
        Logger.getLogger(FastaFileWriter.class.getName()).log(Level.FINER, "writing {0} ({1} bytes)\n", new Object[]{header, Integer.toString(seq.length)});

        BufferedOutputStream out = new BufferedOutputStream(Channels.newOutputStream(channel));
        try {
            out.write('>');
            if (header != null) {
                out.write(header.getBytes("UTF8"));
            }
            out.write('\n');

            if (multiline) {
                for (int i = off, n = i + len; i < n; i += LINE_LENGTH) {
                    out.write(seq, i, Math.min(n - i, LINE_LENGTH));
                    out.write('\n');
                }
            } else {
                for (int i = off, n = i + len; i < n; i += 4096) {
                    out.write(seq, i, Math.min(n - i, 4096));
                }
                out.write('\n');
            }
        } finally {
            out.flush();
        }
        
    }
    
    /**
     * Directly writes the fasta sequence to the file.
     * 
     * @param seq the sequence object to be written
     * 
     * @throws IOException 
     */
    public void write(FastaSequence seq) throws IOException {
        write(seq.header != null ? seq.header : "", seq.multiline, seq.sequence, 0, seq.sequence != null ? seq.sequence.length : 0);        
    }
    
    /**
     * <p>
     * Method to write fasta sequence asynchronously.
     * Method reserves space in the fasta file and 
     * returns a SinkChannel to write the sequence.
     * Method blocks if there is no free thread available in the pool.<br/>
     * <b>NB: Method may fail miserably on Windows when many short sequences written
     * fast - see JDK-6907260. It supposed to be used when sequences are long and 
     * their generation is much slower than disk output speed.</b>
     * </p>
     * 
     * The simple <b>synchronous</b> example:
     * <pre>
     * {@code
     *     FastaSequence seq = ...
     *     try (SinkChannel sink = fastaFileWriter.write(seq)) {
     *         if (seq.sequence != null) {
     *             sink.write(ByteBuffer.wrap(seq.sequence));
     *         }
     *     }
     * }
     * </pre>
     * @param seq
     * @return
     * @throws IOException 
     */
    public SinkChannel write(TFastaSequence seq) throws IOException {
        
        Logger.getLogger(FastaFileWriter.class.getName()).log(Level.FINER, "writing {0} ({1} bytes)\n", new Object[]{seq.header, Integer.toString(seq.length)});
        
        final byte[] header = seq.header != null ? seq.header.getBytes("UTF8") : new byte[] {};
        final ByteBuffer buf = ByteBuffer.allocate(header.length + 2).put((byte)'>').put(header).put((byte)'\n');
        buf.rewind();
        while(channel.write(buf) > 0);

        final long pos = channel.position();
        final long size = seq.multiline ? seq.length + (seq.length / LINE_LENGTH) + 1 : seq.length + 1;

        final ByteBuffer out = channel.map(FileChannel.MapMode.READ_WRITE, pos, size);
        channel.position(pos + size);
                
        final Pipe pipe = Pipe.open(); // JDK-6907260 ?

        if (executor == null) {
            setExecutor(1);
        }
        executor.submit(new FastaSequenceWriter(pipe.source(), out, seq.multiline));

        return pipe.sink();
    }
    
    private void setExecutor(final int threads) {
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

    @Override
    public void close() throws IOException {
        try {
            if (executor != null) {
                executor.shutdown();
                executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
            }
        } catch (InterruptedException ex) {
            throw new IOException(ex.getMessage());
        } finally {
            channel.close();
        }
    }

    public static class FastaSequenceWriter implements Runnable {

        private final SourceChannel src;
        private final ByteBuffer out;
        
        public FastaSequenceWriter(SourceChannel src, ByteBuffer out, boolean multiline) {
            this.src = src;
            this.out = out;
            
            if (multiline) {
                out.limit(LINE_LENGTH);
            }
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
            } finally {
                try {
                    src.close();
                } catch (IOException ex) {
                    Logger.getLogger(FastaFileWriter.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
    }
}
