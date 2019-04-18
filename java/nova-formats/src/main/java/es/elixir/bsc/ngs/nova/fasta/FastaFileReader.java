/**
 * *****************************************************************************
 * Copyright (C) 2017 Spanish National Bioinformatics Institute (INB) and
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

package es.elixir.bsc.ngs.nova.fasta;

import es.elixir.bsc.ngs.nova.gzip.GZipFileInputStream;
import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import static java.nio.file.StandardOpenOption.READ;
import java.util.EnumSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.ZipException;

/**
 * 
 * @author Dmitry Repchevsky
 */

public class FastaFileReader implements Iterable<FastaSequence> {
    
    private final Path path;
    private final boolean gzip;
    private final boolean lazy;

    public FastaFileReader(Path path) throws IOException {
        this(path, false);
    }
    
    /**
     * FastaFileReader constructor.
     * 
     * @param path path to the FASTA file to read.
     * @param lazy if <i>false</i> - do not read the sequence data.
     *             Lazy loading is ignored when FASTA file is gzipped.
     *             Lazy loading may save a memory when one is looking for a
     *             particular read. Do not use it to read all the sequences, especially
     *             when reads are short as it may gravely harm performance.
     * @throws java.io.IOException
     */
    public FastaFileReader(Path path, boolean lazy) throws IOException {
        
        boolean gzip;
        try (GZipFileInputStream in = new GZipFileInputStream(path)) {
            gzip = true;
        } catch (ZipException ex) {
            gzip = false;
        }
        this.gzip = gzip;
        this.path = path;
        this.lazy = lazy;
    }
    
    @Override
    public FastaIterator iterator() {
        return (iterator(null));
    }
    
    public FastaIterator iterator(FastaSequence seq) {
        try {
            InputStream in = gzip ? new GZipFileInputStream(path) :
                    new BufferedInputStream(Files.newInputStream(path, StandardOpenOption.READ));
            return seq == null ? new FastaIterator(in, gzip ? false : lazy) : 
                                 new FastaIterator(in, seq, gzip ? false : lazy);
        } catch (IOException ex) {
            return null;
        }
    }
    
    /**
     * Reads the sequence into the byte buffer
     * 
     * @param buf byte buffer to read the sequence into (may be null)
     * @param seq the sequence to be read into the byte buffer
     * 
     * @return byte buffer where the sequence was put.
     * 
     * @throws java.io.IOException 
     */
    public ByteBuffer read(ByteBuffer buf, FastaSequence seq) throws IOException {
        
        Logger.getLogger(FastaFileReader.class.getSimpleName()).log(Level.FINE, "reading ''{0}'' ({1} bytes)\n", new Object[]{seq.header, seq.length});
        
        if (seq.sequence != null) {
            return buf == null ? ByteBuffer.wrap(seq.sequence) : 
                    buf.put(seq.sequence, 0, Math.min(seq.sequence.length, buf.remaining()));
        }
        
        InputStream in;
        if (gzip) {
            in = new GZipFileInputStream(path);
            in.skip(seq.position);            
        } else {
            FileChannel channel = FileChannel.open(path, EnumSet.of(READ));

            if (seq.multiline) {
                channel.position(seq.position);
                in = new BufferedInputStream(Channels.newInputStream(channel));
            } else {
                try {
                    if (buf == null) {
                        return channel.map(FileChannel.MapMode.READ_ONLY, seq.position, seq.length);
                    } else {
                        final int limit = buf.limit();
                        buf.limit(buf.position() + seq.length);
                        channel.read(buf, seq.position);
                        buf.limit(limit);
                        return buf;
                    }
                } finally {
                    channel.close();
                }
            }
        }

        try {
            if (buf == null) {
                buf = ByteBuffer.allocate(seq.length);
            }

            for (int i = 0, ch = in.read(); i < seq.length; ch = in.read()) {
                if (ch != '\r' && ch != '\n') {
                    i++;
                    buf.put((byte)ch);
                }
            }
        } finally {
            in.close();
        }
        return buf;
    }
}
