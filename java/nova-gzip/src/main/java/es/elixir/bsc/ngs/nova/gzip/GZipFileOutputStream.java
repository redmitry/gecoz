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

package es.elixir.bsc.ngs.nova.gzip;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;

/**
 * @author Dmitry Repchevsky
 */

public class GZipFileOutputStream extends OutputStream 
        implements AutoCloseable {
    
    private int size;
    private long pos;
    private final FileChannel channel;
    private final OutputStream out;
    private final GZipOutputStream gzip;
    
    private final GZipHeader header;

    public GZipFileOutputStream(final Path file) throws IOException {
        
        header = new GZipHeader(0, null, null, 0);
        
        channel = FileChannel.open(file,
                StandardOpenOption.WRITE, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);

        out = new BufferedOutputStream(Channels.newOutputStream(channel));
        gzip = new GZipOutputStream(header, out);
    }

    @Override
    public void write(final int b) throws IOException {
        if (size == 65536) {
            gzip.close();
            
            fix_bsize();
            
            pos = channel.position();

            // write next header;
            header.write(out);
            
            size = 0;
        }
        gzip.write(b);
        size++;
    }

    @Override
    public void write(final byte b[], int off, int len) throws IOException {
        if (len > 0) {
            if (size == 65536) {
                gzip.close();
                fix_bsize();
                pos = channel.position();
                header.write(out);
                size = 0;
            }
            while(true) {
                if (size + len <= 65536) {
                    gzip.write(b, off, len);
                    size += len;
                    break;
                }
                gzip.write(b, off, 65536 - size);
                gzip.close();
                fix_bsize();
                pos = channel.position();
                header.write(out);
                len = len + size - 65536;
                off = off - size + 65536;
                size = 0;
            }                
            
        }
    }

    @Override
    public void close() throws IOException {
        gzip.close();
        fix_bsize();
        channel.close();
    }
    
    private void fix_bsize() throws IOException {
        final short bsize = (short)(channel.position() - pos - 1);
        final ByteBuffer buf = ByteBuffer.allocate(2).order(ByteOrder.LITTLE_ENDIAN);
        buf.putShort(0, bsize);
        while (buf.hasRemaining()) {
            channel.write(buf, buf.position() + pos + 16);
        }
    }
}
