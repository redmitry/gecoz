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

package es.elixir.bsc.ngs.nova.io;

import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * @author Dmitry Repchevsky
 */

public class BitBuffer extends AbstractBitStream {

    public BitBuffer(ByteBuffer buf) {
        super(buf, buf.limit());
    }

    public byte get(long bitIndex) {
        return (byte)(((buf.get((int)(bitIndex >>> 3)) & 0xFF) >> (bitIndex & 7)) & 0x01);
    }
    
    @Override
    public void align() throws IOException {
        buf.position(buf.position() - (bits_left >>> 3));
        bits_left = 0;
    }
}