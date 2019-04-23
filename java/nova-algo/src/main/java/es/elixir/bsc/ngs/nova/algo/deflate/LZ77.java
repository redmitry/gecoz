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

import static es.elixir.bsc.ngs.nova.algo.deflate.Deflater.DEFLATE_WINDOW_SIZE;
import es.elixir.bsc.ngs.nova.algo.string.SAIS;
import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * LZ77 implementation of the DEFLATE Compressed Data Format Specification.
 * The implementation is based on SAIS suffix array algorithm for patterns search.
 * While there is no limit for the data buffer size (apart of the memory),
 * usage of the buffer bigger than 64K slows algorithm as we have to skip patterns 
 * found ahead of compressing point. Encoder is also use heristic to determine best 
 * pattern match to emit (usually closest distance).
 * 
 * @author Dmitry Repchevsky
 */
public class LZ77 {
    
    public final static int MIN_PATTERN_LENGTH = 3;
    public final static int MAX_PATTERN_LENGTH = 258;

    private final LZOutputStream out;

    public LZ77(final LZOutputStream out) {
        this.out = out;
    }

    /**
     * Compress the data. Note that compression is performed from the 
     * current buffer position to the buffer limit.
     * 
     * @param buf the data to be compressed.
     * 
     * @throws IOException 
     */
    public void lz(final ByteBuffer buf) throws IOException {
        
        final int[] sa = new int[buf.limit()];
        final int[] c = SAIS.suffix(buf, sa);
        
        final long[] count = new long[256];
        for (int i = 0, idx = -1; i < c.length; i++) {
            if (c[i] > 0) {
                count[i] = c[i] - idx;
                idx = c[i];
            }
        }
        
        final int[] inv = new int[sa.length];
        final int[] lcp = lcp(buf, sa, inv);

        final DeflateEncodeTable table = new DeflateEncodeTable(count);

        final int est[] = new int[buf.limit() - buf.position() + 1];
        for (int i = 0, j = buf.position(), m = est.length - 1; i < m; i++) {
            est[i + 1] = est[i] + table.bit_lengths[buf.get(i + j) & 0xFF];
        }

        int pos = buf.position();
        do {
            if (lcp[inv[pos]] < MIN_PATTERN_LENGTH && (inv[pos] < lcp.length - 1 && lcp[inv[pos] + 1] < MIN_PATTERN_LENGTH)) {
                out.encode_lit(buf.get(sa[pos]) & 0xFF);
            } else {
                int gain = 0;
                int length = 1;
                int distance = DEFLATE_WINDOW_SIZE; // limit max distance
                
                for (int p = inv[pos] + 1, l = Integer.MAX_VALUE; p < lcp.length && (l = Math.min(l, lcp[p])) >= MIN_PATTERN_LENGTH; p++) {
                    final int d = pos - sa[p];
                    if (d > 0 && d < distance) {
                        final int len_idx = pos - buf.position();
                        final int pattern_bits = est[len_idx + l] - est[len_idx];
                        final int lens_bits = l >= 258 ? 0 : 32 - Integer.numberOfLeadingZeros(l - MIN_PATTERN_LENGTH >>> 3);
                        int total_bits = pattern_bits - lens_bits;
                        if (total_bits < gain) {
                            break;
                        }

                        final int dist_bits = 30 - Integer.numberOfLeadingZeros((d - 1) | 3);

                        total_bits -= dist_bits + 8;
                        if (total_bits > gain) {
                            gain = total_bits;
                            distance = d;
                            length = l;
                        }
                    }
                }

                for (int p = inv[pos] - 1, l = lcp[inv[pos]]; l >= MIN_PATTERN_LENGTH; l = Math.min(l, lcp[p]), p--) {
                    final int d = pos - sa[p];
                    if (d > 0 && d < DEFLATE_WINDOW_SIZE || (l <= length && d < distance)) {
                        final int len_idx = pos - buf.position();
                        final int pattern_bits = est[len_idx + l] - est[len_idx];
                        final int lens_bits = l >= 258 ? 0 : 32 - Integer.numberOfLeadingZeros(l - MIN_PATTERN_LENGTH >>> 3);
                        int total_bits = pattern_bits - lens_bits;
                        if (total_bits < gain) {
                            break;
                        }

                        final int dist_bits = 30 - Integer.numberOfLeadingZeros((d - 1) | 3);

                        total_bits -= dist_bits + 8;
                        if (total_bits > gain) {
                            gain = total_bits;
                            distance = d;
                            length = l;
                        }
                    }
                }

                if (length < MIN_PATTERN_LENGTH) {
                    out.encode_lit(buf.get(sa[pos]) & 0xFF);
                } else {
                    length = Math.min(MAX_PATTERN_LENGTH, distance > 2 && distance < length ? distance : length);

                    out.encode_lens(length--);
                    out.encode_dist(distance);

                    pos += length;
                }
            }
        } while (++pos < sa.length);
    }
    
    private static int[] lcp(final ByteBuffer buf, final int[] sa, final int[] inv) {
        
        for (int i = 0, n = inv.length; i < n; i++) {
            inv[sa[i]] = i;
        }

        int[] lcp = new int[sa.length];
        for (int i = 0, k = 0, n = inv.length; i < n; i++) {
            if (inv[i] == n - 1) {
                k = 0;
                continue;
            }

            int j = sa[inv[i] + 1];
            while (i + k < n && j + k < n &&
                   buf.get(i + k) == buf.get(j + k)) {
                k++;
            }

            lcp[inv[i] + 1] = k;

            if (k > 0) {
                k--;
            }
        }

        return lcp;
    }
}
