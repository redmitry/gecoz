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
import java.io.OutputStream;

/**
 * @author Dmitry Repchevsky
 */

public class DefaultBitOutputStream extends OutputStream implements BitOutputStream {
    
    private long value;
    private int bits_left;

    private final OutputStream out;
    
    public DefaultBitOutputStream(OutputStream out) {
        bits_left = 64;
        this.out = out;
    }
    
    @Override
    public void writeBits(long bits, int nbits) throws IOException {
        if (bits_left > nbits) {
            bits_left -= nbits;
            value = (value << nbits) | bits;
        } else {
            if (bits_left == nbits) {
                bits = bits | (value << nbits);
                bits_left = 64;
                value = 0;
            } else {
                final int tail = nbits - bits_left;
                long tmp = (value << bits_left) | (bits >> tail);
                bits_left = 64 - tail;
                value = bits & (0xFFFFFFFFFFFFFFFFL >>> bits_left);
                bits = tmp;
            }

            for (int i = 56; i >= 0; i -= 8) {
                out.write((byte)(bits >>> i & 0xFF));
            }
        }
    }

    /**
     * Writes the byte to the stream after aligning it
     * 
     * @param b
     * @throws IOException 
     */
    @Override
    public void write(int b) throws IOException {
        flush();
        out.write((byte)(b & 0xFF));
    }
    
    @Override
    public void close() throws IOException {
        flush();
        out.close();
    }
    
    @Override
    public void flush() throws IOException {
        value <<= bits_left;
        for (int i = 56, n = bits_left & 0b1111000; i >= n ; i -= 8) {
            out.write((byte)(value >> i & 0xFF));
        }
        value = 0;
        bits_left = 64;
    }
}
