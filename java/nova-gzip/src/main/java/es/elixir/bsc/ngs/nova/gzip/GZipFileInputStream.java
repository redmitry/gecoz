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

package es.elixir.bsc.ngs.nova.gzip;

import es.elixir.bsc.ngs.nova.algo.deflate.Inflater;
import es.elixir.bsc.ngs.nova.io.DataReaderHelper;
import es.elixir.bsc.ngs.nova.io.FileChannelBitInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.channels.FileChannel;

/**
 * @author Dmitry Repchevsky
 */

public class GZipFileInputStream extends InputStream {
    
    private GZipHeader header;

    private final Inflater inflater;
    private final FileChannelBitInputStream in;
    private final java.util.zip.CRC32 crc;

    private long block_pos;
    
    public GZipFileInputStream(File file) throws IOException {
        final FileChannel channel = new FileInputStream(file).getChannel();
        in = new FileChannelBitInputStream(channel);
        inflater = new Inflater(in, crc = new java.util.zip.CRC32());
        readHeader();
    }

    public long getBlockPosition() {
        return block_pos;
    }
    /**
     * @return the number of uncompressed bytes read from the current gzip chunk.
     * @throws java.io.IOException
     */
    protected int count() throws IOException {
        return inflater.count();
    }

    /**
     * Checks whether the end of stream is reached.
     * 
     * @return -1 if the end of stream is reached, 0 otherwise.
     * 
     * @throws IOException 
     */
    @Override
    public int available() throws IOException {
        if (inflater.finished()) {
            if (header == null || skip(0) < 0) {
                return -1;
            }
        }
        return 0;
    }
    
    /**
     * Moves file position to the concrete position.
     * The position must be the beginning of the gzip chunk
     * 
     * @param pos the position to move
     * @throws IOException 
     */
    protected void move(long pos) throws IOException {
        block_pos = pos;
        in.setPosition(pos);
        inflater.reset();
        readHeader();
    }

    /**
     * Move to the next gzip file block.
     * 
     * @return true if moved next.
     * 
     * @throws java.io.IOException
     */
    public boolean next() throws IOException {
        while (inflater.skip(65536) >= 0);
        return readFooter() >= 0;
    }
    
    /**
     * Skips over and discards n bytes of data from this input stream.
     * 
     * @param n the number of bytes to be skipped.
     * 
     * @return n or -1 if the end of stream is reached.
     * 
     * @throws IOException 
     */
    @Override
    public long skip(final long n) throws IOException {
        long m = n;
        long l;
        while ((l = inflater.skip(m)) < m) {
            if (l < 0) {
                if (readFooter() < 0) {
                    return -1;
                }
            } else {
                m -= l;
            }    
        }
        return n;
    }

    @Override
    public int read() throws IOException {
        final int l = inflater.read();
        if (l < 0 && readFooter() == 0) {
            return read();
        }
        return l;
    }
    
    @Override
    public int read(byte[] buf) throws IOException {
        return read(buf, 0, buf.length);
    }
    
    @Override
    public int read(byte[] buf, int off, int len) throws IOException {
        final int l = inflater.read(buf, off, len);
        if (l < 0 && readFooter() == 0) {
            return read(buf, off, len);
        }
        return l;
    }
    
    private int readFooter() throws IOException {
        long c = DataReaderHelper.readUnsignedInt(in);
        if (c != crc.getValue()) {
            throw new IOException("invalid block crc code");
        }
        crc.reset();
        final long isize = DataReaderHelper.readUnsignedInt(in);
        if (isize == 0) {
            return -1;
        }

        try {
            return readHeader();
        } catch(EOFException ex) {}
        
        return -1; // EOF
    }
    
    private int readHeader() throws IOException {
        block_pos = in.getPosition();
        header = new GZipHeader(in);

        // support for BGZF
        if (header.dsize == 2) {
            header = null;
            return -1;
        }

        inflater.reset();
        return 0;
    }
}
