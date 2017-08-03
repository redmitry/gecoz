package es.elixir.bsc.ngs.nova.algo.ssa;

import es.elixir.bsc.ngs.nova.algo.tree.HuffmanShapedWaveletTree;
import es.elixir.bsc.ngs.nova.algo.tree.IndexWaveletTree;
import es.elixir.bsc.ngs.nova.algo.tree.RankedWTNode;
import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * Sparse Suffix Array Index structure. The index consists in 
 * permuted sparse suffix array stored as Wavelet Tree and a
 * ranked bit vector that marks those BWT symbols for that SA is kept.
 * The marker´s rank (bit '1') corresponds to the position of the sparse index.
 * 
 * @author Dmitry Repchevsky
 */
public class GSSAIndex {
    private final RankedWTNode rank; // bit vector that keeps marked characters (those for which SA is kept)
    private final IndexWaveletTree wsa; // partial suffix array

    public final int sampling_factor; // Sampling Factor 3 =  
    
    /**
     * Create the Sparse Suffix Array Index from its serialized form.
     * The index consists in Permutted Sparse Suffix Array in a form of 
     * Wavelet Matrix and a Ranked Bit Vector located one after another.
     * 
     * @param in
     * @param len
     * @throws IOException 
     */
    public GSSAIndex(ByteBuffer in, long len) throws IOException {

        // having the index size (bytes left in the buffer) and known
        // original length of the Suffix Array, calculate the 
        // sampling factor.
        int _sampling_factor = -1;
        do {
            _sampling_factor++;
        } while (in.remaining() < GSSAIndex.getIndexSize(len, _sampling_factor));
            
        sampling_factor = _sampling_factor;
        
        rank = new RankedWTNode(in, len);
        wsa = new IndexWaveletTree(in, (len + (1 << sampling_factor) - 1) >> sampling_factor);
    }
    
    public GSSAIndex(HuffmanShapedWaveletTree tree, int sampling_factor) throws IOException {
        this(tree, null, sampling_factor);
    }

    /**
     * Creates the Sparse Suffix Array Index from a BWT Wavelet Tree.
     * Because of almost random BWT access, the time is similar to the de-novo
     * SSA construction.
     * 
     * @param tree
     * @param c
     * @param sampling_factor
     * 
     * @throws IOException 
     */
    public GSSAIndex(HuffmanShapedWaveletTree tree, long[] c, int sampling_factor) throws IOException {
        this.sampling_factor = sampling_factor;

        if (c == null) {
            c = new long[256];
            for (int i = 0; i < 256; i++) {
                c[i] = tree.occ(i, tree.shape.length); // NB!!! will it fail on chars that are not in the tree?
            }
        }

        final int[] sa = new int[(int)tree.shape.length]; // 
        
        final long mask = 0xFFFFFFFFFFFFFFFFL >>> (64 - sampling_factor);
        
//        long idx = -1;
//        long i = tree.length - 1;
//        while (i >= 0) {
//            final long rs = tree.getRS(idx <= bwt ? idx + 1 : idx);
//            idx = (int)(c[(int)rs] + (rs >>> 32));
//            if ((i-- & mask) == 0) {
//                sa[(int)idx] = i + 1;
//            }
//        }

        int[] ssa = new int[(sa.length + (1 << sampling_factor) - 1) >> sampling_factor];
        
        rank = new RankedWTNode(sa.length);
//        for (int i = 0, j = 0, n = sa.length; i < n; i++) {
//            int pos = sa[i];
//            if (pos-- > 0) {
//                ssa[j++] = pos >> sampling_factor;
//                rank.put(1);
//            } else {
//                rank.put(0);
//            }
//        }
        rank.flush();
        
        wsa = new IndexWaveletTree(ssa);
    }
    
    private GSSAIndex(int[] sa, int sampling_rate, ByteBuffer out) throws IOException {
    
        sampling_factor = 31 - Integer.numberOfLeadingZeros(sampling_rate);
        
        int[] ssa = new int[(sa.length + (1 << sampling_factor) - 1) >> sampling_factor];
        
        final int mask = 0xFFFFFFFF >>> (32 - sampling_factor);

        rank = new RankedWTNode(out, sa.length);
        for (int i = 0, j = 0, n = sa.length; i < n; i++) {
            final int pos = sa[i];
            if ((pos & mask) == 0) {
                ssa[j++] = pos >> sampling_factor;
                rank.put(1);
            } else {
                rank.put(0);
            }
        }
        rank.flush();

        wsa = new IndexWaveletTree(ssa, out);
    }

    public long size() {
        return rank.size;
    }
    
    /**
     * Returns the Suffix Array index for the SSA position.
     * Because this Suffix Array is sparse, some positions have no index stored.
     * 
     * @param pos the position to get Suffix Array index
     * 
     * @return the Suffix Array index or Integer.MIN if no index stored at the position
     */
    public long get(long pos) {
        return rank.get(pos) == 0 ? Integer.MIN_VALUE : wsa.get(rank.count(pos) - 1) << sampling_factor;
    }

    /**
     * Finds the position for the given Suffix Array index.
     * Because this Suffix Array is sparse, some indexes are absent.
     * 
     * @param idx the suffix array index
     * 
     * @return the position for the index or Integer.MIN if no index stored
     */
    public long find(long idx) {
        final long sidx = idx >> sampling_factor;
        return idx == sidx << sampling_factor ? rank.findOne(wsa.find(sidx) + 1) : Integer.MIN_VALUE;
    }
    
    public static GSSAIndex write(int[] sa, int sampling_rate, ByteBuffer out) throws IOException {
        return new GSSAIndex(sa, sampling_rate, out);
    }
    /**
     * Calculates the index size for a Suffix Array.
     * 
     * @param size the size of the indexed Suffix Array
     * @param sampling_factor the sampling factor which is a power of 2 (1,2,3,4 ...)
     * 
     * @return the size of the index (Sparse Suffix Array + Ranks Bit Vector)
     */
    public static long getIndexSize(long size, int sampling_factor) {
        final long ssa_len = (size + (1 << sampling_factor) - 1) >> sampling_factor;
        final long ssa_size = IndexWaveletTree.size(ssa_len);
        final long rnk_size = RankedWTNode.bytes(size);
        return ssa_size + rnk_size;
    }
}
