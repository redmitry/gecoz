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

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;


/**
 * The BitInputStream implementation based on memory mapped file.
 * 
 * @author Dmitry Repchevsky
 */

public class FileChannelBitInputStream extends InputStream implements BitInputStream {
    
    private final FileChannel channel;
    private ByteBuffer buf;

    private int total_bits;
    private byte bits_left;
    private long value;

    public FileChannelBitInputStream(FileChannel channel) throws IOException {
        this(channel, 64 * 1024 * 1024);
    }
    
    public FileChannelBitInputStream(FileChannel channel, int size) throws IOException {
       this.channel = channel;
       map(size);
    }
    
    private void map(int size) throws IOException {
        final long remain = channel.size() - channel.position();
        if (remain > 0) {
            int len = (int)Math.min(remain, size);
            if (len >= 8) {
                len = len - len % 8;
                total_bits = len * 8;

                buf = channel.map(FileChannel.MapMode.READ_ONLY, channel.position(), len);
                buf.order(ByteOrder.LITTLE_ENDIAN);
            } else {
                total_bits = len * 8;

                buf = ByteBuffer.allocate(8);
                buf.order(ByteOrder.LITTLE_ENDIAN);
                channel.read(buf);
                buf.rewind();
            }
        }
    }

    public long getPosition() throws IOException {
        return channel.position() + this.buf.position() - (bits_left >>> 3);
    }

    public void setPosition(long pos) throws IOException {
        bits_left = 0;
        channel.position(pos);
        map(buf.limit());
    }
    
    @Override
    public final long peekBits(final int nbits) throws IOException {
        if (bits_left == 0) {
            if (total_bits < nbits) {
                channel.position(channel.position() + buf.limit());
                map(buf.limit());
            }

            bits_left = 64;
            value = buf.getLong();
            return value;
        }
        
        final int last = (int)(value >>> (64 - bits_left));
        if (bits_left >= nbits) {
            return last;
        }

        final int next;
        if (total_bits >= nbits) {
            next = buf.getShort(buf.position());
        } else {
            final ByteBuffer buffer = ByteBuffer.allocate(2);
            buffer.order(ByteOrder.LITTLE_ENDIAN);
            
            while (buffer.hasRemaining() && 
                   channel.read(buffer, channel.position() + buf.limit() + buffer.position()) >= 0) {
            }
            next = buffer.getShort(0);
        }
        return last | (next << bits_left);
    }
    
    @Override
    public final void skipBits(final int nbits) throws IOException { 
        bits_left -= nbits;
        if (bits_left < 0) {
            if (total_bits < nbits) {
                channel.position(channel.position() + buf.limit());
                map(buf.limit());
                total_bits -= bits_left;
            }
            bits_left += 64;
            value = buf.getLong();
        }
        total_bits -= nbits;
    }
    
    @Override
    public final long readBits(final int nbits) throws IOException {
        if (bits_left == 0) {
            if (total_bits < nbits) {
                channel.position(channel.position() + buf.limit());
                map(buf.limit());
            }
            total_bits -= nbits;
            value = buf.getLong();
            bits_left = (byte)(64 - nbits);
            return value;
        } else {
            long bits = value >>> (64 - bits_left);
            if (bits_left < nbits) {
                if (total_bits < nbits) {
                    channel.position(channel.position() + buf.limit());
                    map(buf.limit());
                    total_bits += bits_left;
                }
                value = buf.getLong();
                bits |= (value << bits_left);
                bits_left += (64 - nbits);
            } else {
                bits_left -= nbits;
            }
            total_bits -= nbits;
            return bits;
        }
    }
    
    @Override
    public final long skip(final long nbytes) throws IOException {
        align();

        this.buf.position(this.buf.position() - (bits_left >>> 3));
        bits_left = 0;
        
        final int remaining = this.buf.remaining();
        if (nbytes <= remaining) {
            total_bits -= nbytes * 8;
            buf.position(this.buf.position() + (int)nbytes);
            return nbytes;
        }
        
        final long position = channel.position();
        
        channel.position(position + this.buf.position() + nbytes);
        map(this.buf.limit());
        
        return channel.position() - position + this.buf.limit() - remaining; // CHECK!!!
    }

    @Override
    public final void align() throws IOException {
        skipBits(bits_left & 0x7); // align bytes
    }
    
    @Override
    public final int read() throws IOException {
        align();
        
        if (total_bits >= 8) {
            return (int) (readBits(8) & 0xFF);
        }
        
        if (channel.size() == channel.position()) {
            return -1; // EOF
        }

        // assert bits_left == 0
        channel.position(channel.position() + buf.limit());
        map(buf.limit());
        
        return read();
    }
    
    
    @Override
    public final int read(byte[] buf, final int off, final int len) throws IOException {
        align();
        
        this.buf.position(this.buf.position() - (bits_left >>> 3));
        bits_left = 0;
        
        final int remaining = this.buf.remaining();
        if (len <= remaining) {
            this.buf.get(buf, off, len);
            total_bits -= len * 8;
            return len;
        } 

        if (remaining == 0 && 
            channel.position() >= channel.size()) {
            return -1;
        }
        
        this.buf.get(buf, off, remaining);
        
        total_bits -= remaining * 8;
        
        channel.position(channel.position() + this.buf.limit());
        map(this.buf.limit());

        return remaining;
    }
}
