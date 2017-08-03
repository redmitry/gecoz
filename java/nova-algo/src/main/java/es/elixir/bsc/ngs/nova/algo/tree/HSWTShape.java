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

import es.elixir.bsc.ngs.nova.algo.deflate.DeflateEncodeTable;
import es.elixir.bsc.ngs.nova.algo.deflate.DeflateLengthsTable;
import es.elixir.bsc.ngs.nova.algo.deflate.DeflateLookupTable;
import es.elixir.bsc.ngs.nova.io.BitBuffer;
import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * @author Dmitry Repchevsky
 */

public class HSWTShape {
    
    /**
     * The size in bytes of the RFC 1951 3.2.7 Deflate Length Table.
     */
    public final long size;
    
    /**
     * The length of the data (sum of all counts[])
     */
    public final long length;
    public final long[] counts;

    public final DeflateEncodeTable encode;
    public final DeflateLookupTable decode;
    
    public HSWTShape(long[] counts) {
        
        this.counts = counts;
        
        encode = new DeflateEncodeTable(counts);
        decode = new DeflateLookupTable(encode.bit_lengths);

        int[] lengths = new int[256]; // the size of bit verctors (in bits)
        
        long len = 0;
        for (int i = 0; i < 256; i++) {
            if (counts[i] > 0) {
                len += counts[i];
                final int code = encode.table[i];
                for (int j = 0, n = encode.bit_lengths[i]; j < n; j++) {
                    int idx = code & (0x0000FFFF >>> (16 - j));
                    idx |= (0x8000 >>> (15 - j));
                    idx = decode.getSymbol(idx);
                    lengths[idx] += counts[i];
                }
            }
        }
        
        long sz = (DeflateLengthsTable.length(encode.bit_lengths) + 7) >>> 3;
        for (int i = 0; i < 256; i++) {
            if (lengths[i] > 0) {
                sz += RankedWTNode.bytes(lengths[i]);
            }
        }
        
        length = len;
        size = sz;
    }
    
    private HSWTShape(ByteBuffer in, long length) throws IOException {
        
        this.length = length;
        
        final int pos = in.position();
        
        BitBuffer buf = new BitBuffer(in);
        DeflateLengthsTable lengths = new DeflateLengthsTable(buf, 256);        
        buf.align();
        
        size = in.position() - pos;

        encode = new DeflateEncodeTable(lengths.d_tree);
        decode = new DeflateLookupTable(lengths.d_tree);
        
        counts = null;
    }
    
    public static HSWTShape read(ByteBuffer in, long length) throws IOException {
        return new HSWTShape(in, length);
    }

    public void write(ByteBuffer out) throws IOException {
        BitBuffer buf = new BitBuffer(out);
        new DeflateLengthsTable(encode.bit_lengths).write(buf);
        buf.flush();
    }
}
