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

package es.elixir.bsc.ngs.nova.bam;

import java.io.File;
import java.io.IOException;
import java.util.zip.DataFormatException;
import es.elixir.bsc.ngs.nova.gzip.GZipFileInputStream;

/**
 * @author Dmitry Repchevsky
 */

public class BAMFileInputStream extends GZipFileInputStream {

    private final BAMHeader header;
    
    public BAMFileInputStream(File fbam) throws IOException, DataFormatException {
        this(fbam, null);
    }
    
    public BAMFileInputStream(File fbam, File fbai) throws IOException, DataFormatException {
        super(fbam);
        header = new BAMHeader(this);
    }
    
    public int getRefCount() {
        return header.refs.length;
    }
    
    /**
     * Gets the current position index
     * 
     * @return
     * @throws IOException 
     */
    public long index() throws IOException {
        return count() | (super.getBlockPosition() << 16);
    }
    
    /**
     * Moves the position to the BAM index.
     * NB: The method changes the contract of the parent accepting the compound BAM index 
     * and not just a file offset.
     * 
     * @param index
     * @throws IOException 
     */
    @Override
    public void move(long index) throws IOException {
        super.move(index >>> 16);
        index &= 0xFFFF;
        for (long i = index, j = 0; i > 0 && (j = skip(i)) >= 0; i -= j);
    }
}
