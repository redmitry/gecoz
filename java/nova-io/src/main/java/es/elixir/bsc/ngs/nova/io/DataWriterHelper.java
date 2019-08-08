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

package es.elixir.bsc.ngs.nova.io;

import java.io.IOException;
import java.io.OutputStream;

/**
 * @author Dmitry Repchevsky
 */

public class DataWriterHelper {
    
    public static void writeUnsignedShort(final OutputStream out, final int value) throws IOException {
        out.write(value & 0xFF);
        out.write((value >> 8) & 0xFF);
    }
    
    public static void writeUnsignedInt(final OutputStream out, final long value) throws IOException {
        out.write((int)(value & 0xFF));
        out.write((int)((value >> 8) & 0xFF));
        out.write((int)((value >> 16) & 0xFF));
        out.write((int)((value >> 24) & 0xFF));
    }
    
    public static void writeLong(final OutputStream out, final long value) throws IOException {
        writeUnsignedInt(out, value & 0xFFFFFFFFL);
        writeUnsignedInt(out, value >>> 32);
    }
}
