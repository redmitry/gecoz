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

import java.io.EOFException;
import java.io.Flushable;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * @author Dmitry Repchevsky
 */

public abstract class AbstractBitStream implements BitInputStream, BitOutputStream, Flushable {
    
    public final long size;
    
    protected long value;
    protected byte bits_left;
    protected ByteBuffer buf;
    
    public AbstractBitStream(ByteBuffer buf, long size) {
        this.buf = buf;
        this.buf.order(ByteOrder.LITTLE_ENDIAN);
        
        this.size = size;
    }

    public void put(int bit) throws IOException {
        writeBits(bit, 1);
    }

    public void write(ByteBuffer out) throws IOException {
        flush();
        final int position = buf.position();
        buf.position(0);
        out.put(buf);
        buf.position(position);
    }
    
    @Override
    public long peekBits(int nbits) throws IOException {

        if (bits_left == 0) {
            value = getLong();
            bits_left = 64;
            return value;
        }
        
        final long last = value >>> (64 - bits_left);
        if (bits_left >= nbits) {
            return last;
        }

        final long next = getLong(buf.position());
        return last | (next << bits_left);
    }

    @Override
    public void skipBits(int nbits) throws IOException {
        bits_left -= nbits;
        if (bits_left < 0) {
            if (buf.hasRemaining()) {
                bits_left += 64;
                value = getLong();
            } else {
                throw new EOFException();
            }
        }
    }

    @Override
    public long readBits(int nbits) throws IOException {

        if (bits_left == 0) {
            value = getLong();
            bits_left = (byte)(64 - nbits);
            return value;
        } else {
            long bits = value >>> (64 - bits_left);
            if (bits_left < nbits) {
                value = getLong();
                bits |= (value << bits_left);
                bits_left += (64 - nbits);
            } else {
                bits_left -= nbits;
            }
            return bits;
        }
    }

    public void rewind() {
        if (bits_left > 0) {
            putLong(value);
            bits_left = 0;
        }
        buf.rewind();
    }

    @Override
    public void writeBits(long bits, int nbits) throws IOException {

        if (bits_left > nbits) {
            value |= bits << (64 - bits_left);
            bits_left -= nbits;
        } else if (bits_left == 0) {
            value = bits;
            bits_left = (byte) (64 - nbits);
        } else if (bits_left < nbits) {
            putLong(value | (bits << (64 - bits_left)));
            value = bits >> bits_left;
            bits_left += (64 - nbits);
        } else {
            putLong(value | bits << (64 - bits_left));
            bits_left = 0;
        } 
    }

    /**
     * Method flushes buffered long value to the internal ByteBuffer.
     * It results to the byte alignment.
     * 
     * @throws IOException 
     */
    @Override
    public void flush() throws IOException {
        if (bits_left > 0) {
            final int pos = buf.position();
            final int len = (71 - bits_left) >>> 3;
            putLong(value);
            buf.position(pos + len);
            bits_left = 0;
        }
    }

    protected long getLong() {
        if (buf.remaining() < Long.BYTES) {
            // can't read entire Long, so read byte by byte
            long l = buf.get() & 0xFF;
            for (int i = 8; buf.hasRemaining(); i += 8) {
                l |= (buf.get() & 0xFFL) << i;
            }
            return l;
        }
        return buf.getLong();
    }
    
    protected long getLong(long index) {
        if (buf.limit() - index < Long.BYTES) {
            // can't read entire Long, so read byte by byte
            long l = 0;
            for (int i = 0, n = buf.limit(); index < n; i += 8, index++) {
                l |= (buf.get((int)index) & 0xFFL) << i;
            }
            return l;
        }
        return buf.getLong((int)index);
    }
    
    protected void putLong(long value) {
        if (buf.remaining() < Long.BYTES) {
            // can't write entire Long, so write what we can
            for (int i = 0, n = 64 - bits_left; buf.hasRemaining() && i <= n; i+=8) {
                buf.put((byte)((value >>> i) & 0xFF));
            }
        } else {
            buf.putLong(value);
        }
    }
}
