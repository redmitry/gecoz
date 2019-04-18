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

import es.elixir.bsc.ngs.nova.io.BitOutputStream;
import java.io.IOException;

/**
 * Deflater Block Header class used for the compression.
 * 
 * @author Dmitry Repchevsky
 */

public class DeflaterBlockHeader extends BlockHeader {

    public final DeflaterDynHeader header;
    
    public DeflaterBlockHeader(final boolean bfinal) {
        super(HUFF_FIXED, bfinal);
        this.header = null;
    }

    public DeflaterBlockHeader(final int len, final boolean bfinal) {
        super(NO_COMPRESSION, bfinal);
        this.len = len;
        this.header = null;
    }
    
    public DeflaterBlockHeader(final DeflaterDynHeader header,
                               final boolean bfinal) {
        
        super(HUFF_DYNAMIC, bfinal);
        this.header = header;
    }

    public final void write(final BitOutputStream out) throws IOException {
        out.writeBits((btype << 1) | (bfinal ? 1 : 0), 3);
        switch(btype) {
            case NO_COMPRESSION: out.writeBits(0, 5); // align
                                 out.writeBits(len & 0xFFFF, 16);
                                 out.writeBits(len ^ 0xFFFF, 16);
            case HUFF_FIXED:     break;
            case HUFF_DYNAMIC:   header.write(out);
                                 break;
            default: throw new IOException("invalid btype");
        }
    }
}
