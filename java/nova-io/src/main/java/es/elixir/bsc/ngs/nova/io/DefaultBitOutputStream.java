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

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * @author Dmitry Repchevsky
 */

public class DefaultBitOutputStream extends FilterOutputStream implements BitOutputStream {
    
    private long value;
    private byte bits_pushed;

    public DefaultBitOutputStream(final OutputStream out) {
        super(out);
    }

    @Override
    public void writeBits(final long bits, final int nbits) throws IOException {
        value |= bits << bits_pushed;
        bits_pushed += nbits;
        if (bits_pushed >= 64) {
            put(value);
            bits_pushed -= 64;
            value = bits >>> (nbits - bits_pushed);
        }
    }

    @Override
    public void align() throws IOException {
        if (bits_pushed > 0) {
            for (int i = 0; i < bits_pushed; i += 8) {
                write((int)((value >>> i) & 0xFF));
            }

            bits_pushed = 0;
            value = 0;
        }        
    }

    @Override
    public void flush() throws IOException {
        align();
        out.flush();
    }
    
    @Override
    public void close() throws IOException {
        flush();
    }

    private void put(final long value) throws IOException {
        for (int i = 0; i < 64; i += 8) {
            write((int)((value >>> i) & 0xFF));
        }
    }
}
