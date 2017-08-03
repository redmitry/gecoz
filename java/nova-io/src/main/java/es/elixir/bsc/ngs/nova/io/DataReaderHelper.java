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

package es.elixir.bsc.ngs.nova.io;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

/**
 * @author Dmitry Repchevsky
 */

public class DataReaderHelper {
    
    public static int readUnsignedShort(InputStream in) throws IOException {
        final int b0 = in.read();
        final int b1 = in.read();
        if (b1 < 0) {
            throw new EOFException();
        }
        return b0 | (b1 << 8);
    }
    
    public static long readUnsignedInt(InputStream in) throws IOException {
        final int b0 = in.read();
        final int b1 = in.read();
        final int b2 = in.read();
        final int b3 = in.read();
        if (b3 < 0) {
            throw new EOFException();
        }
        return (b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)) & 0xFFFFFFFFL;
    }
    
    public static long readLong(InputStream in) throws IOException {
        final long l0 = readUnsignedInt(in);
        final long l1 = readUnsignedInt(in);
        
        return l0 | (l1 << 32);
    }
}
