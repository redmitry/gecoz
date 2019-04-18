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

package es.elixir.bsc.ngs.nova.algo.deflate;

import static es.elixir.bsc.ngs.nova.algo.deflate.Deflater.DEFLATE_WINDOW_SIZE;
import es.elixir.bsc.ngs.nova.io.BitOutputStream;
import java.io.Flushable;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;

/**
 * RFC1951 Deflater.
 * 
 * @author Dmitry Repchevsky
 */

public class DeflaterOutputStream extends OutputStream implements Flushable, AutoCloseable {
    
    private int size;
    private final byte[] data;
    private final ByteBuffer buf;

    private final LZ77 lz;
    private final Deflater deflater;
    
    public DeflaterOutputStream(final BitOutputStream out) {
        data = new byte[DEFLATE_WINDOW_SIZE * 2];
        buf = ByteBuffer.wrap(data);
        deflater = new Deflater(out);
        lz = new LZ77(deflater);
    }
    
    @Override
    public void write(final int value) throws IOException {
        if (size == data.length) {
            buf.limit(data.length);
            lz.lz(buf);
            System.arraycopy(data, data.length - DEFLATE_WINDOW_SIZE, data, 0, DEFLATE_WINDOW_SIZE);
            size = DEFLATE_WINDOW_SIZE;
            buf.position(DEFLATE_WINDOW_SIZE);
        }
        
        data[size++] = (byte)value;
    }

    @Override
    public void write(final byte b[], int off, int len) throws IOException {
        if (size + len > data.length) {
            while (len > data.length - size) {
                System.arraycopy(b, off, data, size, data.length - size);
                buf.limit(data.length);
                lz.lz(buf);
                System.arraycopy(data, data.length - DEFLATE_WINDOW_SIZE, data, 0, DEFLATE_WINDOW_SIZE);
                len -= data.length - size;
                off += data.length - size;
                size = DEFLATE_WINDOW_SIZE;
                buf.position(DEFLATE_WINDOW_SIZE);
            }
        }
        System.arraycopy(b, off, data, size, len);
        size += len;
    }

    /**
     * Flushes current data into the new deflate block.
     * Each block stores its own Huffman statistics.
     * LZ search is also performed over previous block.
     * 
     * @throws IOException 
     */
    @Override
    public void flush() throws IOException {
        deflate(false);
        if (size > DEFLATE_WINDOW_SIZE) {
            System.arraycopy(data, size - DEFLATE_WINDOW_SIZE, data, 0, DEFLATE_WINDOW_SIZE);
            size = DEFLATE_WINDOW_SIZE;
        }
        buf.position(size);
    }
    
    /**
     * Flushes current data into the new deflate block and closes deflate stream
     * (writing 'last block' marker).
     * Method does not close underlying stream.
     * 
     * @throws IOException 
     */
    @Override
    public void close() throws IOException {
        deflate(true);
        size = 0;
        buf.position(0);
    }
    
    private void deflate(final boolean bfinal) throws IOException {
        if (size > buf.position()) {
            buf.limit(size);
            lz.lz(buf);
        }
        deflater.deflate(bfinal);
    }
}
