/**
 * *****************************************************************************
 * Copyright (C) 2015 Spanish National Bioinformatics Institute (INB) and
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

package es.elixir.bsc.ngs.nova.algo.huffman;

import es.elixir.bsc.ngs.nova.io.BitOutputStream;
import java.io.IOException;

/**
 * The class generates Huffman codes - lower frequency node is given '1'.
 * The codes are inverted (lower to higher bits = root to leaves)
 * 
 * @author Dmitry Repchevsky
 */

public class HuffmanEncodeTable {

    public final short[] table;
    public final byte[] bit_lengths;
    
    protected HuffmanEncodeTable(final byte[] bit_lengths) {
        this.bit_lengths = bit_lengths;
        table = new short[bit_lengths.length];
    }

    public HuffmanEncodeTable(final long[] counts) {

        table = new short[counts.length];
        bit_lengths = new byte[counts.length];
        int[] bt = new int[counts.length]; // tree circular pointers

        for (int i = 1, n = counts.length; i < n; i++) {
            int idx1 = 0;
            int idx2 = 0;
            long min1 = Long.MAX_VALUE, min2 = Long.MAX_VALUE;
            for (int j = 0; j < n; j++) {
                final long fq = counts[j];
                if (fq > 0) {
                    if (fq < min1) {
                        idx2 = idx1;
                        min2 = min1;
                        idx1 = j;
                        min1 = fq;
                    } else if (fq < min2) {
                        idx2 = j;
                        min2 = fq;
                    }
                }
            }

            if (min1 == Long.MAX_VALUE || min2 == Long.MAX_VALUE) {
                break;
            }

            counts[idx1] = Long.MIN_VALUE;
            counts[idx2] = min1 + min2;
            
            if (bt[idx1] == 0) {
                bt[idx1] = bt[idx2] < 0 ? bt[idx2] : ~idx2;
                bt[idx2] = ~idx1;
            } else if (bt[idx2] == 0) {
                bt[idx2] = bt[idx1];
                bt[idx1] = ~idx2;
            } else {
                // merge two branches
                final int idx = bt[idx1];
                bt[idx1] = bt[idx2];
                bt[idx2] = idx;
            }
            
            // increase all bits counters in the circle;
            int idx = idx1;
            do {
                idx = ~bt[idx]; 
                table[idx] <<= 1;
                bit_lengths[idx]++;
            } while (idx != idx2);
            do {
                idx = ~bt[idx];
                table[idx] = (short)((table[idx] << 1) | 1);
                bit_lengths[idx]++;
            } while (idx != idx1);
        }
    }

    public final void putSymbol(int symbol, final BitOutputStream out) throws IOException {
        out.writeBits(table[symbol], bit_lengths[symbol]);
    }
}
