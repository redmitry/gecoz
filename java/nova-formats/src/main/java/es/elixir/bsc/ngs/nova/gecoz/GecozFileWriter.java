/**
 * *****************************************************************************
 * Copyright (C) 2017 Spanish National Bioinformatics Institute (INB) and
 * Barcelona Supercomputing Center
 *
 * Modifications to the initial code base are copyright of their respective
 * authors, or their employers as appropriate.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 *****************************************************************************
 */

package es.elixir.bsc.ngs.nova.gecoz;

import es.elixir.bsc.ngs.nova.algo.ssa.GSSAIndex;
import es.elixir.bsc.ngs.nova.algo.string.SuffixArray;
import es.elixir.bsc.ngs.nova.algo.tree.HSWTShape;
import es.elixir.bsc.ngs.nova.algo.tree.HuffmanShapedWaveletTree;
import es.elixir.bsc.ngs.nova.algo.tree.HuffmanShapedWaveletTree.DataSource;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import static java.nio.file.StandardOpenOption.CREATE;
import static java.nio.file.StandardOpenOption.READ;
import static java.nio.file.StandardOpenOption.WRITE;
import static java.nio.file.StandardOpenOption.TRUNCATE_EXISTING;
import java.util.EnumSet;
import java.io.Closeable;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.RunnableFuture;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Dmitry Repchevsky
 */

public class GecozFileWriter implements Closeable {
    
    private ExecutorService executor;
    
    private final FileChannel ref_channel;
    private final FileChannel ssa_channel;
    
    private final int sampling_rate;
    
    public GecozFileWriter(Path path) throws IOException {
        this(path, null);
    }
    
    public GecozFileWriter(Path path, int sampling_rate) throws IOException {
        this(path, null, sampling_rate);
    }

    public GecozFileWriter(Path ref_path, Path ssa_path) throws IOException {
        this(ref_path, ssa_path, 32);
    }
    
    public GecozFileWriter(Path ref_path, Path ssa_path, int sampling_rate) throws IOException {
        this(ref_path, ssa_path, sampling_rate, 1);
    }
    
    /**
     * 
     * @param ref_path
     * @param ssa_path
     * @param sampling_rate the sampling rate for the SSA index (8,16,32...(
     * @param th the desired number of threads to use
     * @throws IOException 
     */
    public GecozFileWriter(Path ref_path, Path ssa_path, int sampling_rate, int th) throws IOException {
        ref_channel = FileChannel.open(ref_path, EnumSet.of(CREATE,READ,WRITE, TRUNCATE_EXISTING));
        
        if (ssa_path == null) {
            PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:**.gcz");
            String ssa_fname = ref_path.getFileName().toString();
            if (matcher.matches(ref_path)) {
                ssa_fname = ssa_fname.substring(0, ssa_fname.length() - 3);
            }
            ssa_path = ref_path.resolveSibling(ssa_fname + "gcx");
        }
        
        ssa_channel = FileChannel.open(ssa_path, EnumSet.of(CREATE,READ,WRITE, TRUNCATE_EXISTING));
        
        this.sampling_rate = sampling_rate;
        
        final int threads = Math.min(th,Runtime.getRuntime().availableProcessors());
        
        Logger.getLogger(GecozFileWriter.class.getName()).log(Level.FINER, "writer uses {0} threads\n", threads);

        executor = new WriterPoolExecutor(threads);
    }
    
    /**
     * Adds a new block with a generalized string
     * 
     * @param headers
     * @param in
     * @throws IOException 
     */
    public void write(String[] headers, ByteBuffer in) throws IOException {

        // calculate characters' frequencies
        final long counts[] = new long[256];
        for (int i = 0, n = in.limit(); i < n; i++) {
            counts[in.get(i) & 0xFF]++;
        }
        
        HSWTShape shape = new HSWTShape(counts);
        
        // total block size
        long ref_pos = ref_channel.position();
        long ref_block_size = GecozRefBlockHeader.getBlockHeaderLength(headers) + shape.size;

        ByteBuffer out = ref_channel.map(FileChannel.MapMode.READ_WRITE, ref_pos, ref_block_size);
        out.order(ByteOrder.LITTLE_ENDIAN);
        
        GecozRefBlockHeader ref_header = new GecozRefBlockHeader(headers, ref_block_size, in.remaining());
        ref_header.write(out);

        ref_channel.position(ref_pos + ref_block_size);

        final long idx_pos = ssa_channel.position();
        final long idx_size = GSSAIndex.getIndexSize(in.remaining(), 31 - Integer.numberOfLeadingZeros(sampling_rate));
        final long idx_block_size = GecozSSABlockHeader.getBlockHeaderLength() + idx_size;
        
        ByteBuffer idx = ssa_channel.map(FileChannel.MapMode.READ_WRITE, idx_pos, idx_block_size);
        idx.order(ByteOrder.LITTLE_ENDIAN);
        
        GecozSSABlockHeader ssa_header = new GecozSSABlockHeader(headers, idx_size);
        ssa_header.write(idx);

        ssa_channel.position(idx_pos + idx_block_size);

        executor.submit(new BlockWriter(in, out, idx, shape, sampling_rate));
    }

    @Override
    public void close() throws IOException {
        try {
            executor.shutdown();
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch (InterruptedException ex) {
            throw new IOException(ex.getMessage());
        } finally {
            ref_channel.close();
            ssa_channel.close();
        }
    }
    
    public static class WriterPoolExecutor extends ThreadPoolExecutor
                                           implements RejectedExecutionHandler {

        public WriterPoolExecutor(final int threads) {
            super(1, threads, 0L, TimeUnit.MILLISECONDS, new ArrayBlockingQueue<>(1, true));
            setRejectedExecutionHandler(this);
        }

        @Override
        protected <T> RunnableFuture<T> newTaskFor(Runnable runnable, T value) {
            return new BlockWriterFutureTask(runnable);
        }
        
        @Override
        public final void rejectedExecution(Runnable runnable, ThreadPoolExecutor executor) {
            if (!executor.isShutdown()) {
                /* The queue is full (the executor has rejected the task).
                 * Put the task in the queue directly blocking the thread
                 * which submitted the task (called write() method).
                 */
                try {
                    getQueue().put(runnable);
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    Logger.getLogger(GecozFileWriter.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
            
        @Override
        public void afterExecute(Runnable runnable, Throwable th) {
            super.afterExecute(runnable, th);

            final BlockWriterFutureTask task = (BlockWriterFutureTask)runnable;
            try {
                task.get();
            } catch(Exception ex) {
                final int sz = getMaximumPoolSize() - 1;
                if (sz == 0) {
                    // can't process even in one thread...
                    ex.printStackTrace(System.err);
                    System.exit(-1);
                }

                try {
                    getQueue().put(new BlockWriterFutureTask(task.runnable));
                } catch (InterruptedException ex2) {
                    Thread.currentThread().interrupt();
                    Logger.getLogger(GecozFileWriter.class.getName()).log(Level.SEVERE, null, ex2);
                }
                setMaximumPoolSize(sz);
            }
        }
    }

    /**
     * Custom FutureTask which provides an access to the Runnable.
     */
    public static class BlockWriterFutureTask extends FutureTask {
        final Runnable runnable;
        public BlockWriterFutureTask(Runnable runnable) {
            super(runnable, null);
            this.runnable = runnable;
        }
    }

    public static class BlockWriter implements Runnable {

        private final ByteBuffer in;
        private final ByteBuffer out;
        private final ByteBuffer idx;
        private final HSWTShape shape;
        private final int sampling_rate;
        
        public BlockWriter(ByteBuffer in, ByteBuffer out, ByteBuffer idx, HSWTShape shape, int sampling_rate) {
            this.in = in;
            this.out = out;
            this.idx = idx;
            this.shape = shape;
            this.sampling_rate = sampling_rate;
        }

        @Override
        public void run() {
            Logger.getLogger(GecozFileWriter.class.getName()).log(Level.FINE, "indexing {0} bytes\n", in.limit());
            
            try {
                final int[] sa = new int[in.limit()];
                SuffixArray.suffix(in, sa);
                
                ExecutorService exs = Executors.newSingleThreadExecutor();
                exs.submit(() -> {
                    try {
                        shape.write(out);
                        HuffmanShapedWaveletTree.write(shape, new BWTDataSource(in, sa), out);
                    } catch(IOException ex) {
                        throw new RuntimeException(ex);
                    }
                });

                GSSAIndex.write(sa, sampling_rate, idx);

                exs.shutdown();
                exs.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
            } catch(OutOfMemoryError ex) {
                Logger.getLogger(GecozFileWriter.class.getName()).log(Level.WARNING, "warning: low memory (free: {0} bytes)\n", Runtime.getRuntime().freeMemory());
                throw new RuntimeException(ex);
            } catch(IOException | InterruptedException ex) {
                throw new RuntimeException(ex);
            }
        }
    }
    
    /**
     * HuffmanShapedWaveletTree.DataSource
     */
    public static class BWTDataSource implements DataSource {

        private final ByteBuffer in;
        private final int[] sa;
        
        public BWTDataSource(final ByteBuffer in, final int[] sa) {
            this.in = in;
            this.sa = sa;
        }
        
        @Override
        public final byte get(int idx) {
            return sa[idx] == 0 ? in.get(in.limit() - 1) : in.get(sa[idx] - 1);
        }

        @Override
        public final int length() {
            return in.limit();
        }
    }
}
