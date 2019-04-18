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

import es.elixir.bsc.ngs.nova.algo.deflate.DeflaterOutputStream;
import es.elixir.bsc.ngs.nova.io.DefaultBitOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * @author Dmitry Repchevsky
 */

public class GZipOutputStream extends DeflaterOutputStream {

    private long size;
    private final java.util.zip.CRC32 crc;
    private final OutputStream out;

    public GZipOutputStream(final OutputStream out) throws IOException {
        this(new GZipHeader(), out);
    }
    
    public GZipOutputStream(final GZipHeader header, 
                            final OutputStream out) throws IOException {
        super(new DefaultBitOutputStream(out));
        
        this.out = out;
        crc = new java.util.zip.CRC32();
        header.write(out);
    }
    
    private GZipOutputStream(final GZipHeader header, 
                             final OutputStream out, 
                             final DefaultBitOutputStream bits) throws IOException {
        super(bits);
        
        this.out = bits;
        crc = new java.util.zip.CRC32();
        header.write(out);
    }

    @Override
    public void write(final int b) throws IOException {
        super.write(b);
        crc.update(b);
        size++;
    }

    @Override
    public void write(final byte b[]) throws IOException {
        this.write(b, 0, b.length);
    }

    @Override
    public void write(final byte b[], int off, int len) throws IOException {
        super.write(b, off, len);
        crc.update(b, off, len);
        size += len;
    }
    
    /**
     * Closes gzip stream writing the last deflate block and 
     * gzip footer (CRC32 + ISIZE) without closing underlying output stream.
     * 
     * @throws IOException 
     */
    @Override
    public void close() throws IOException {
        super.close();

        // write crc
        final long crc32 = crc.getValue();
        out.write((int)(crc32 & 0xFF));
        out.write((int)((crc32 >>> 8) & 0xFF));
        out.write((int)((crc32 >>> 16) & 0xFF));
        out.write((int)((crc32 >>> 24) & 0xFF));
        
        out.write((int)(size & 0xFF));
        out.write((int)((size >>> 8) & 0xFF));
        out.write((int)((size >>> 16) & 0xFF));
        out.write((int)((size >>> 24) & 0xFF));
        
        out.flush();
        
        size = 0;
        crc.reset();
    }
}
