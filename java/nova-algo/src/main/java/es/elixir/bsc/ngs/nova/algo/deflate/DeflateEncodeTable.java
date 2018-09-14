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

package es.elixir.bsc.ngs.nova.algo.deflate;

import es.elixir.bsc.ngs.nova.algo.huffman.HuffmanEncodeTable;
import java.util.Arrays;

/**
 * The class generates 'deflate' Huffman codes for the alphabet.
 * By default it generates 'deflate' codes - all the codes of the same lengths
 * are lexically ordered. Note that 'deflate' codes are reversed.
 * 
 * @author Dmitry Repchevsky
 */

public class DeflateEncodeTable extends HuffmanEncodeTable {
    
    private final static int MAX_BITS = 15;
    
    public DeflateEncodeTable(final long[] counts) {
        super(Arrays.copyOf(counts, counts.length));
        restrict_lengths(counts);
        remap_codes();
    }
    
    public DeflateEncodeTable(final byte[] bit_lengths) {
        super(bit_lengths);
        remap_codes();
    }
    
    private void restrict_lengths(final long[] counts) {

        long count = 0;

        final int[] bl_count = new int[bit_lengths.length];
        for (int i = 0, n = bit_lengths.length; i < n; i++) {
            final int bl = bit_lengths[i];
            if (bl > 0) {
                count += bl;
                bl_count[bl]++;
            }
        }

        if (count <= 1) {
            return; // special case where all characters are the same
        }

        int nodes = 1; // non-leaf nodes at the level 1;
        for (int i = 1; i <= MAX_BITS && nodes > 0; i++) {
            nodes <<= 1;
            if (bl_count[i] != 0) {
                nodes -= bl_count[i];
            }
        }

        // if no non-leaf nodes left - all bit lengths <= MAX_BITS
        if (nodes > 0) {
            // calculate how many nodes we have to reallocate
            
            nodes = -nodes; // we already have some 'empty' non-leaf nodes
            for (int i = 0, n = bit_lengths.length; i < n; i++) {
                if (bit_lengths[i] > MAX_BITS) {
                    bit_lengths[i] = MAX_BITS;
                    nodes++;
                }
            }
            
            // sort by bit_lengths / counts and keep the index in the last (lowest) byte
            
            // 0xLL_CCCCCCCC_IIII (I - index, C - count, L - bit length)
            long[] list = new long[counts.length];
            for (int i = 0, n = list.length; i < n; i++) {
                list[i] = ((long)bit_lengths[i] << 48) | ((long)counts[i] << 16) | (long)i;
            }

            Arrays.sort(list);

            do {
                loop:
                for (int i = MAX_BITS - 1; i > 0; i--) {
                    for (int level = i; level < MAX_BITS; level++) {
                        long bit_length = (long)(level + 1) << 48;
                        for (int j = 0, n = list.length; j < n; j++) {
                            final byte bl = (byte)(list[j] >>> 48);
                            if (bl == level) {
                                list[j] = list[j] & 0xFF00FFFFFFFFFFFFL | bit_length;
                                
                                // moving the leaf top-down releases a node on the leaf level
                                // what gives us 2^n leaves at the MAX_BITS level
                                nodes -= 1 << (MAX_BITS - 1 - level);
                                if (nodes <= 0) {
                                    break loop;
                                }
                            }
                        }
                    }
                }
                
                for (int level = MAX_BITS; nodes < 0 && level > 0; level--) {
                    long bit_length = (long)(level - 1) << 48;
                    for (int i = list.length - 1; nodes < 0 && i >= 0; i--) {
                        final byte bl = (byte)(list[i] >>> 48);
                        if (bl == level) {
                            list[i] = list[i] & 0xFF00FFFFFFFFFFFFL | bit_length;
                            nodes += 1 << (MAX_BITS - level);
                        }
                    }
                }
            } while (nodes != 0);
            
            // put bit_lengths back
            for (int i = 0, n = bit_lengths.length; i < n; i++) {
                bit_lengths[(int)(list[i] & 0xFF)] = (byte)(list[i] >>> 48);
            }
        }
    }
    
    private void remap_codes() {
        final int[] bl_count = new int[MAX_BITS + 1];
        for (int i = 0, n = bit_lengths.length; i < n; i++) {
            final int bl = bit_lengths[i];
            if (bl > 0) {
                bl_count[bl]++;
            }
        }

        final int[] next_code = new int[MAX_BITS + 1];

        for (int bits = 1, code = 0; bits <= MAX_BITS; bits++) {
            code = (code + bl_count[bits - 1]) << 1;
            next_code[bits] = code;
        }

        for (int i = 0; i < bit_lengths.length; i++) {
            int len = bit_lengths[i];
            if (len != 0) {
                table[i] = (short)(reverse(next_code[len]) >> (16-len));
                next_code[len]++;
            }
        }        
    }
    
    private int reverse(int i) {
        i = (i & 0x00005555) << 1 | (i >>> 1) & 0x00005555;
        i = (i & 0x00003333) << 2 | (i >>> 2) & 0x00003333;
        i = (i & 0x000000F0F) << 4 | (i >>> 4) & 0x00000F0F;
        return (i >>> 8) | (i << 8) & 0xFFFF;
    }
}
