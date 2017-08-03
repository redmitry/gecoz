/**
 * *****************************************************************************
 * Copyright (C) 2017 Spanish National Bioinformatics Institute (INB) and
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

package es.elixir.bsc.ngs.nova.gecoz;

import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.zip.DataFormatException;

/**
 *
 * @author Dmitry Repchevsky
 */

public class GecozSSABlockHeader {
    public final static String MAGIC = "GecozSSA";
    
    public final byte version = 1;
    public final long len;
    public final long hash;
    
    public GecozSSABlockHeader(String[] headers, long len) {
        this.len = len;
        this.hash = GecozRefBlockHeader.getBlockHeaderHash(headers);
    }
    
    public GecozSSABlockHeader(ByteBuffer buf) throws IOException, DataFormatException {
        if (buf.remaining() < 25) {
            throw new EOFException();
        }

        if (buf.getLong() != 0x4153537A6F636547L |
            buf.get() != version) { // "ASSzoceG" - LITTLE ENDIAN
            throw new DataFormatException();
        }
        
        
        this.len = buf.getLong();
        this.hash = buf.getLong();
    }
    
    public int getBlockLength() {
        return getBlockHeaderLength();
    }
    
    public void write(ByteBuffer buf) {
        buf.put(MAGIC.getBytes()); // 8 bytes
        buf.put(version);          // 1 byte
        buf.putLong(len);          // 8 bytes
        buf.putLong(hash);         // 8 bytes
    }
    
    public static int getBlockHeaderLength() {
        return 25;
    }
}
