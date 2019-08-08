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

package es.elixir.bsc.ngs.nova.bam;

import es.elixir.bsc.ngs.nova.gzip.GZipFileOutputStream;
import java.io.IOException;
import java.nio.file.Path;

/**
 * <p>
 * 
 * </p>
 * 
 * @author Dmitry Repchevsky
 */

public class BAMFileOutputStream extends GZipFileOutputStream {
    
    /**
     * <p>
     * Create a new BAM file and write the BAM header.
     * </p>
     * 
     * @param file   file path for the new BAM file
     * @param header BAM header for the BAM file
     * 
     * @throws IOException 
     */
    public BAMFileOutputStream(final Path file, final BAMHeader header) throws IOException {
        super(file);
        
        header.write(this);
        flush(); // flush b
    }

    public void write(final BAMRecord record) throws IOException {
        record.write(this);
    }
}
