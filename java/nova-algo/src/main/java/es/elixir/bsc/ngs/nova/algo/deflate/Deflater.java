/**
 * *****************************************************************************
 * Copyright (C) 2019 Spanish National Bioinformatics Institute (INB) and
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

import static es.elixir.bsc.ngs.nova.algo.deflate.Inflater.BLC_DIST;
import static es.elixir.bsc.ngs.nova.algo.deflate.Inflater.BLC_DIST_EXT_BITS;
import static es.elixir.bsc.ngs.nova.algo.deflate.Inflater.BLC_LENS;
import static es.elixir.bsc.ngs.nova.algo.deflate.Inflater.BLC_LENS_EXT_BITS;
import es.elixir.bsc.ngs.nova.io.BitOutputStream;
import java.io.IOException;
import java.util.Arrays;

/**
 * Huffman encoding implementation of the DEFLATE Compressed Data Format Specification.
 * Deflater accumulates LZ77 data and then flushes it into the BitOutputStream.
 * All data saved as one deflate block. Deflater may be reused multiple times.
 * While compressing, deflater checks whether there is a gain storing a pointer and 
 * duplicate the string if appropriate.
 * 
 * @author Dmitry Repchevsky
 */

public class Deflater implements LZOutputStream {
    
    public final static int DEFLATE_WINDOW_SIZE = 32768;
        
    private final BitOutputStream out;

    private final InflaterOutput inflater;
    
    private int sz_lit;
    private int sz_ext;
    private int sz_dist;
    
    private short[] lit;
    private short[] ext;
    private byte[] dist;
    
    private final long[] count_lit;
    private final long[] count_dist;
    
    public Deflater(final BitOutputStream out) {
        this.out = out;
        
        inflater = new InflaterOutput(null);
         
        count_lit = new long[286];
        count_dist = new long[BLC_DIST.length];
        
        lit = new short[DEFLATE_WINDOW_SIZE];
        ext = new short[DEFLATE_WINDOW_SIZE];
        dist = new byte[DEFLATE_WINDOW_SIZE];
    }

    @Override
    public void encode_lit(final int value) throws IOException {
        if (sz_lit == lit.length) {
            lit = Arrays.copyOf(lit, sz_lit << 1);
        }
        lit[sz_lit++] = (short)value;
        
        count_lit[value]++;
    }

    @Override
    public void encode_lens(final int value) throws IOException {
        
        final int lens_bits;
        final int lens_code;
        if (value == 258) {
            lens_bits = 0;
            lens_code = 28;
        } else {
            lens_bits = 32 - Integer.numberOfLeadingZeros(value - 3 >>> 3);
            lens_code = (lens_bits << 2) + ((value - 3 >>> lens_bits) & 7);
        }

        encode_lit(lens_code + 257);

        if (lens_bits > 0) {
            if (sz_ext == ext.length) {
                ext = Arrays.copyOf(ext, sz_ext << 1);
            }
            ext[sz_ext++] = (short)(value - BLC_LENS[lens_code]);
        }      
    }

    @Override
    public void encode_dist(final int value) throws IOException {
        if (sz_dist == dist.length) {
            dist = Arrays.copyOf(dist, sz_dist << 1);
        }

        final int dist_bits = 30 - Integer.numberOfLeadingZeros((value - 1) | 3);
        final int dist_code = (dist_bits << 1) + ((value - 1 >>> dist_bits) & 3);

        dist[sz_dist++] = (byte)dist_code;
        count_dist[dist_code]++;

        if (dist_bits > 0) {
            if (sz_ext == ext.length) {
                ext = Arrays.copyOf(ext, sz_ext << 1);
            }
            ext[sz_ext++] = (short)(value - BLC_DIST[dist_code]);
        }
    }

    /**
     * Writes stored data into the deflate block along with the block header
     * 
     * @param bfinal - whether this block is the last one in deflate stream.
     * 
     * @throws IOException 
     */
    public void deflate(final boolean bfinal) throws IOException {

        count_lit[256]++;
                
        final DeflaterDynHeader header = new DeflaterDynHeader(
                new DeflateEncodeTable(count_lit), new DeflateEncodeTable(count_dist));
        new DeflaterBlockHeader(header, bfinal).write(out);

        for (int i = 0, j = 0, k = 0; i < sz_lit; i++) {
            final int symbol = lit[i];
            
            if (symbol > 256) {
                final int len_code = symbol - 257;
                final int ext_len_bits = BLC_LENS_EXT_BITS[len_code];
                final int ext_len = ext_len_bits == 0 ? 0 : ext[j++];
                final int len = BLC_LENS[len_code] + ext_len;
                final int dist_code = dist[k++];
                final int ext_dist_bits = BLC_DIST_EXT_BITS[dist_code];
                final int ext_dist = ext_dist_bits == 0 ? 0 : ext[j++];
                final int dis = BLC_DIST[dist_code] + ext_dist;
                
                /*** check if there is no gain on the pattern compression ****/
                int gain = 0;
                for (int d = dis, n = dis - len; d > n; d--) {
                    final byte sym_len = header.litLenTable.bit_lengths[inflater.peek(d)];
                    if (sym_len == 0) {
                        // we can't expand the character which is not in the huffman table.
                        // the copied pattern should be in the previous deflate block.
                        gain = Integer.MAX_VALUE;
                        break;
                    }
                    gain += sym_len;
                }

                gain -= header.litLenTable.bit_lengths[symbol];
                gain -= ext_len_bits;
                gain -= header.distLenTable.bit_lengths[dist_code];
                gain -= ext_dist_bits;

                if (gain < 0) {
                    for (int l = 0; l < len; l++) {
                        final int sym = inflater.peek(dis);
                        inflater.push(sym);
                        header.litLenTable.putSymbol(sym, out);
                    }
                    continue;
                }
                /*************************************************************/
                
                inflater.inflate(dis, len);
                
                header.litLenTable.putSymbol(symbol, out);
                if (ext_len_bits > 0) {
                    out.writeBits(ext_len, ext_len_bits);
                }

                header.distLenTable.putSymbol(dist_code, out);
                if (ext_dist_bits > 0) {
                    out.writeBits(ext_dist, ext_dist_bits);
                }
            } else {
                header.litLenTable.putSymbol(symbol, out);
                inflater.push(symbol);
            }
        }

        header.litLenTable.putSymbol(256, out);

        if (bfinal) {
            out.align();
            inflater.reset();
        }

        sz_lit = 0;
        sz_ext = 0;
        sz_dist = 0;
        
        Arrays.fill(count_lit, 0);
        Arrays.fill(count_dist, 0);
    }
}
