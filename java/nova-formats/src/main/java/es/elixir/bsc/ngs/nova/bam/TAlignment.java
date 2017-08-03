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

import java.io.IOException;
import java.io.InputStream;
import java.util.zip.DataFormatException;
import es.elixir.bsc.ngs.nova.io.DataReaderHelper;

/**
 * Truncated BAM Alignment object that reflects only the alignment block size, 
 * reference sequence ID and the alignment position.
 * 
 * @author Dmitry Repchevsky
 */

public class TAlignment {
    protected int block_size;
    protected int refID;
    protected int pos;

    public TAlignment(InputStream in) throws IOException, DataFormatException {
        block_size = (int)DataReaderHelper.readUnsignedInt(in);
        refID = (int)DataReaderHelper.readUnsignedInt(in);
        pos = (int)DataReaderHelper.readUnsignedInt(in) + 1;
    }
}