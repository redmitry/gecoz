package es.elixir.bsc.ngs.nova.algo.ssa;

import es.elixir.bsc.ngs.nova.algo.tree.HuffmanShapedWaveletTree;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;

/**
 * Generalized Succinct Suffix Array implementation
 * 
 * @author Dmitry Repchevsky
 */

public class GSSA {
    
    private long[] c; // ends of buckets' positions in SA
    private long[] e; // sorted ends of the strings ('\0' positions)
    
    private GSSAIndex index;
    private final HuffmanShapedWaveletTree tree; // wavelet tree that keeps the bwt

    /**
     * Create Generalized Succinct Suffix Array from existing tree and index.
     * 
     * @param tree
     * @param index
     * @throws IOException 
     */
    public GSSA(HuffmanShapedWaveletTree tree, GSSAIndex index) throws IOException {
        this.tree = tree;
        this.index = index;
    }
        
    public long getLength() {
        return tree.shape.length;
    }

    /**
     * Get the length of a string.
     * 
     * @param nstr the index of a string to get the length from.
     * 
     * @return the length of nstr string.
     */
    public long getLength(int nstr) throws IOException {
        index();
        
        if (nstr < 0 || nstr >= e.length) {
            throw new IndexOutOfBoundsException("String index " + nstr + " is out of bound");
        }
        
        return nstr == 0 ? e[nstr] : e[nstr] - e[nstr - 1] - 1;
    }
    
    /**
     * Extracts the original string S from the Succinct Suffix Array.
     * 
     * @param buf the buffer where the extracted string is put.
     * @param nstr the nth string to extract from.
     * @param from the position of the nth string to start the extraction.
     * 
     * @throws IOException 
     */
    public void extract(ByteBuffer buf, int nstr, long from) throws IOException {
        index();

        if (nstr < 0 || nstr >= e.length) {
            throw new IndexOutOfBoundsException("String index " + nstr + " is out of bound");
        }
        
        if (nstr > 0) {
            // from index for the nth string
            from += e[nstr - 1] + 1;
        }

        // the last character position in the generalized string 'S'
        // from which we are going to perform the backward extract
        long pos = Math.min(e[nstr], from + buf.remaining()) - 1;

        // the position we have SA index kept
        final long sapos = ((pos >> index.sampling_factor) + 1) << index.sampling_factor;

        // the BWT position 
        long idx = sapos < tree.shape.length ? index.find(sapos) : 0;

        // skip off characters
        long n = Math.min(sapos, tree.shape.length - 1) - pos;
        while (--n > 0) {
            long rs = tree.getRS(idx);
            idx = (int)(c[(int)rs] + (rs >>> 32));            
        }

        final int bpos = buf.position() + (int)(pos - from);
        for (int i = bpos, m = buf.position(); i >= m; i--) {
            long rs = tree.getRS(idx);
            buf.put(i, (byte)(rs & 0xFF));
            idx = (int)(c[(int)rs] + (rs >>> 32));            
        }
        buf.position(bpos + 1);
    }

    /**
     * Count the number of the string occurrences in the SSA.
     * 
     * @param str the string to find
     * @return the number of string occurrences found
     * 
     * @throws IOException 
     */
    public long[] count(byte[] str) throws IOException {
        long[][] sa = find(str);
        
        if (sa == null || sa.length == 0) {
            return null;
        }

        long[] count = new long[sa.length];
        for (int i = 0, n = sa.length; i < n; i++) {
            count[i] = sa[i] != null && sa[i].length > 0 ? sa[i].length : 0;
        }
        return count;
    }

    /**
     * Finds the string in the SSA.
     * 
     * @param str the string to find
     * @return the array[nth][pos] with all occurrences of the string.
     *         it returns null if no matches found at all and nulls 
     *         for any nth string with no matches (array[nth] == null)
     * 
     * @throws IOException 
     */
    public long[][] find(byte[] str) throws IOException {
        long[] sa = search(str);
        
        if (sa == null || sa.length == 0) {
            return null;
        }

        // find the correspondence between 'S' position matches and substrings
        // it takes log(n)log(m) where 'n' is a number of matches (sa.length) and
        // 'm' is a number of substrings (e.length)
        Arrays.sort(sa);
        
        long[][] res = new long[e.length][];
        for (int i = 0, idx1 = 0, n = e.length; i < n; i++) {
            final int idx2 = -Arrays.binarySearch(sa, idx1, sa.length, e[i]) - 1;
            if (idx2 > idx1) {
                final long pos = i > 0 ? e[i-1] + 1 : 0;
                res[i] = new long[idx2 - idx1];
                for (int j = 0, m = res[i].length; j < m; j++) {
                    res[i][j] = sa[idx1 + j] - pos;
                }     
                idx1 = idx2;
            }
        }
        return res;
    }

    private long[] search(byte[] str) throws IOException {
        index();

        int ch = (byte)(str[str.length - 1] & 0xFF);
        long sp = c[ch];
        long ep = ch < c.length - 1 ? c[ch + 1] - 1 : tree.shape.length - 1;
        for (int i = str.length - 2; sp <= ep && i >= 0; i--) {
            ch = (byte)(str[i] & 0xFF);
            sp = c[ch] + tree.occ(ch, sp - 1) + 1;
            ep = c[ch] + tree.occ(ch, ep);
        }

        if (ep < sp) {
            return new long[0]; // not found
        }
        
        long[] res = new long[(int)(ep - sp + 1)];
        for (int i = 0, n = res.length; i < n; i++) {
            res[i] = locate(sp++);
        }
        return res;
    }

    /**
     * Reconstructs SSA index from the BWT.
     * 
     * @throws IOException 
     */
    private void index() throws IOException {
        if (c == null) {
            long idx = tree.shape.length;
            c = new long[256];
            for (int i = 255; i >= 0; i--) {
                final long rank = tree.occ(i, tree.shape.length - 1);
                if (rank >= 0) {
                    idx -= rank + 1;
                }
                c[i] = idx;
            }
        }
        
        if (index == null) {
            index = new GSSAIndex(tree, c, 4);
        }
        
        if (e == null) {
            e = new long[(int)c[1]];
            for (int i = 0, n = e.length; i < n; i++) {
                e[i] = locate(i);
            }
            Arrays.sort(e);
        }
    }
    
    private long locate(long idx) {
        long len = 0;
        long sa = index.get(idx);
        while(sa < 0) {
            len++;
            final long rs = tree.getRS(idx);
            idx = (int)(c[(int)rs] + (rs >>> 32));
            sa = index.get(idx);
        }
        return sa + len;
    }
}
