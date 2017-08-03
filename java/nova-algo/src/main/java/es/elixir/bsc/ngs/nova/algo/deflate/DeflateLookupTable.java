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

import es.elixir.bsc.ngs.nova.io.BitInputStream;
import java.io.IOException;

/**
 * @author Dmitry Repchevsky
 */

public class DeflateLookupTable {
    
    private final static int MAX_BITS = 15;
    private final short[] table;

    public DeflateLookupTable(final byte[] bit_lengths) {
        
        final int[] bl_count = new int[MAX_BITS + 1];
        for (int i = 0, n = bit_lengths.length; i < n; i++) {
            final int bl = bit_lengths[i];
            if (bl > 0) {
                bl_count[bl]++;
            }
        }

        final int[] next_code = new int[MAX_BITS + 1];

        int tree_size = 512;
        for (int bits = 1, code = 0; bits <= MAX_BITS; bits++) {
            final int count = bl_count[bits];
            final int tail = code << 25 >>> (42 - bits);
            code += count << (16 - bits);
            if (bits > 9) {
                // Nº permutations in main table (only first 9 bits counts)
                // multiply "extra" bits ^ 2
                tree_size += ((code - next_code[bits - 1]) << (bits - 9) >>> 7) + tail;
            }
            next_code[bits] = code;
        }
        
        table = new short[tree_size];

        if (tree_size > 512) {
            
            // put pointers to the extended tables
            
            for (int bits = 10, ptr = 512, step = 64; bits <= MAX_BITS; bits++, step >>>= 1) {
                final int ext = bits - 9;
                final int code = next_code[bits - 1];
                final int next = next_code[bits];
                final int tail = code << 25 >>> (42 - bits);
                if (tail > 0) {
                    ptr+= tail; // reserve space for smaller bit lengths permutations
                    
                    // fix number of extension bits to read (first 4 bits)
                    final int idx = reverse(code) & 0b111111111;
                    table[idx] = (short)(table[idx] & 0xFFFFFFF0 | ext);
                }
                
                for (int i = code; i < next; i += step, ptr++) {
                    final int idx = reverse(i) & 0b111111111;
                    if (table[idx] == 0) {
                        table[idx] = (short)(-ptr << 4 | ext);
                    }
                }
            }
        }
        
        for (int i = 0; i < bit_lengths.length; i++) {
            final int bits = bit_lengths[i];
            if (bits > 0) {
                final int code = next_code[bits - 1];
                int revcode = reverse(code);
                if (bits < 9) {
                    for (int j = revcode, l = 0x1 << bits; j < 512; j += l) {
                        table[j] = (short) (i << 4 | bits);
                    }
                } else if (bits == 9) {
                    table[revcode] = (short) (i << 4 | bits);
                } else {
                    final short ptr = table[revcode & 0x1FF];
                    final int ext = ptr & 0b1111;
                    int idx = (code << 25 >>> (32 - ext)) - (ptr >> 4);  // subtable index
                    for (int k = idx + (1 << (ext + 9 - bits)); idx < k; idx++) {
                        table[idx] = (short)(i << 4 | (bits - 9));
                    }
                }
                next_code[bits - 1] = code + (0x01 << (16 - bits));
            }
        }
    }
    
    /**
     * Decodes a symbol from the bit stream.
     * 
     * @param in - deflate bit stream
     * @return the decoded symbol
     * @throws IOException 
     */
    public final int getSymbol(final BitInputStream in) throws IOException {
        final int peek = (int)(in.peekBits(9) & 0b111111111);
        int symbol = table[peek];
        final int bits = symbol & 0b1111;
        if (symbol >= 0) {
            in.skipBits(bits); // skip actually read bits;
        } else {
            in.skipBits(9);
            final int idx = (reverse((int)in.peekBits(bits)) >>> (16 - bits)) - (symbol >> 4);
            symbol = table[idx];
            in.skipBits(symbol & 0b1111);
        }
        return symbol >>> 4;
    }

    /**
     * Decodes a symbol from the Huffman code.
     * 
     * @param code - the Huffman code
     * @return the decoded symbol
     */
    public final int getSymbol(final int code) {
        final int peek = (int)(code & 0b111111111);
        int symbol = table[peek];
        if (symbol < 0) {
            final int idx = ((reverse(code >>> 9)) >>> (16 - (symbol & 0b1111))) - (symbol >> 4);
            symbol = table[idx];
        }
        return symbol >>> 4;
    }

    /**
     * Decodes a symbol from the Huffman code.
     * 
     * @param code - the Huffman code
     * @param nbits - the number of bits in the code
     * @return the decoded symbol or Integer.MIN_VALUE if not enough bits 
     */
    public final int getSymbol(final int code, final int nbits) {
        final int peek = (int)(code & 0b111111111);
        int symbol = table[peek];
        int len = symbol & 0b1111;
        if (symbol < 0) {
            len += 9;
            final int idx = ((reverse(code >>> 9)) >>> (16 - len)) - (symbol >> 4);
            symbol = table[idx];
        }
        
        return nbits >= len ? symbol >>> 4 : Integer.MIN_VALUE;
    }

    private int reverse(int i) {
        i = (i & 0x00005555) << 1 | (i >>> 1) & 0x00005555;
        i = (i & 0x00003333) << 2 | (i >>> 2) & 0x00003333;
        i = (i & 0x000000F0F) << 4 | (i >>> 4) & 0x00000F0F;
        return (i >>> 8) | (i << 8) & 0xFFFF;
    }
}
