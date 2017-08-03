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

import es.elixir.bsc.ngs.nova.io.BitBuffer;
import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * @author Dmitry Repchevsky
 */

public class DefaultWTNode extends BitBuffer implements WaveletTreeNode {
 
    public DefaultWTNode(ByteBuffer buf) {
        super(buf);
    }
    
    /**
     * Creates the new Wavelet Tree Node
     * 
     * @param size - the size of the node in bits;
     */
    public DefaultWTNode(final int size) {
        super(ByteBuffer.allocate((size + 7) >>> 3));
    }

    @Override
    public long size() {
        return size;
    }

    @Override
    public void put(int bit) throws IOException {
        writeBits(bit, 1);
    }

    @Override
    public void write(ByteBuffer out) throws IOException {
        flush();
        final int position = buf.position();
        buf.position(0);
        out.put(buf);
        buf.position(position);
    }

    @Override
    public long count(long idx) {
        long count = 0;

        long n =  idx >>> 6;
        for (int i = 0; i < n; i++) {
            count += Long.bitCount(getLong(i << 3));
        }

        return count + Long.bitCount(getLong(n << 3) << (63 - (idx & 63)));
    }
}
