/**
 * *****************************************************************************
 * Copyright (C) 2016 Spanish National Bioinformatics Institute (INB) and
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

package es.elixir.bsc.ngs.nova.algo.tree;

import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * <p>
 * Wavelet Tree that represents an index array. The index array must be a 
 * permutation of numbers from 0 to n, where n is the size of the array.
 * Operations take log(n) time.
 * </p>
 * 
 * @author Dmitry Repchevsky
 */

public class IndexWaveletTree {
    
    private final RankedWTNode[] nodes;
    
    /**
     * Create IndexWaveletTree from the index array
     * 
     * @param sa
     * @throws IOException 
     */
    public IndexWaveletTree(int[] sa) throws IOException {
        this(sa, null);
    }
    
    /**
     * <p>Create IndexWaveletTree from its serialized form.</p>
     * <p>
     * Technically, the serialized form is just a chain of <i>log(n) + 1</i> sized
     * binary vectors where <i>n</i> is the size of the index.
     * For example an array of 64 indexes is represented as log(64)+1 bit vectors
     * of length 64 (occupies 7 longs).
     * <p>
     * 
     * @param in buffer that contains chained binary vectors.
     * @param size size of the index.
     */
    public IndexWaveletTree(ByteBuffer in, long size) {
        int hibit = 64 - Long.numberOfLeadingZeros(size);
        nodes = new RankedWTNode[hibit];
        
        while (hibit-- > 0) {
            nodes[hibit] = new RankedWTNode(in, size);
        }
    }
    
    /**
     * Create IndexWaveletTree from the index array using ByteBuffer as a storage
     * 
     * @param sa
     * @param dst
     * @throws IOException 
     */
    public IndexWaveletTree(int[] sa, ByteBuffer dst) throws IOException {
        int[] _ssa = new int[sa.length];
        int hibit = 32 - Integer.numberOfLeadingZeros(sa.length);
        nodes = new RankedWTNode[hibit];//(ssa.length * hibit);
        
        while (hibit-- > 0) {
            nodes[hibit] = dst == null ? new RankedWTNode(sa.length) : new RankedWTNode(dst, sa.length);
            final int mask = 0xFFFFFFFF << hibit;
            for (int i = 0, n = sa.length; i < n; i++) {
                final int idx = sa[i];
                final int block = idx & mask; // block start
                int c = Math.min(block + (1 << hibit), sa.length) - 1; // where to store counter
                int ptr = _ssa[c];
                if (ptr >= 0) {
                    _ssa[c] = ~block;
                    _ssa[block] = idx;
                } else {
                    _ssa[c] = --ptr;
                    _ssa[~ptr] = idx;
                }
                byte bit = (byte)((idx >> hibit) & 1);
                nodes[hibit].put(bit);
            }

            nodes[hibit].flush();
            int[] tmp = sa;
            sa = _ssa;
            _ssa = tmp;
        }
    }
    
    /**
     * Get the size of the index.
     * @return the size of the index.
     */
    public long getSize() {
        return nodes[0].size;
    }
    
    /**
     * Get index stored at the position.
     * @param pos position of the index.
     * @return the index found at the position {@code pos}.
     */
    public long get(long pos) {        
        long code = 0;
        
        for (int i = 63 - Long.numberOfLeadingZeros(getSize()), block = 0; i >= 0; i--) {
            final int bit = nodes[i].get(pos);
            long bits = nodes[i].count(pos);
            code = code << 1 | bit;
            if (bit == 0) {
                bits = pos - bits - (block >>> 1);
            } else {
                bits -= (block >>> 1) + 1;
                block += (1 << i);
            }
            pos = block + bits;
        }
        
        return code;
    }
    
    /**
     * Find the position of the index.
     * 
     * @param idx the index to search.
     * @return the position of the index.
     */
    public long find(long idx) {
        long pos = 0;
        
        for (int i = 0, n = 64 - Long.numberOfLeadingZeros(getSize()); i < n; i++) {
            final long bit = (idx >>> i) & 1;
            final long block = idx & (0xFFFFFFFFFFFFFFFEL << i); // block start
            pos = bit == 0 ? nodes[i].findZero((block >>> 1) + pos + 1, block, Math.min(block + (2 << i), nodes[i].size) - 1) :
                             nodes[i].findOne((block >>> 1) + pos + 1, block, Math.min(block + (2 << i), nodes[i].size) - 1);
            
            pos -= block;
        }
        
        return pos;
    }
    
    /**
     * Calculates the space in bytes needed to store <i>len</i> indexes.
     * 
     * @param len the size of index array (number of indexes).
     * @return the size in bytes needed to store the index array of length <i>len</i>.
     */
    public static long size(long len) {
        return RankedWTNode.bytes(len) * (64L - Long.numberOfLeadingZeros(len));
    }
}
