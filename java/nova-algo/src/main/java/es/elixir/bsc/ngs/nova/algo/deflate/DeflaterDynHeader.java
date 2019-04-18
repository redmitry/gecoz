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
 * @author Dmitry Repchevsky
 */

public class DeflaterDynHeader {
    
    public final DeflateEncodeTable litLenTable;
    public final DeflateEncodeTable distLenTable;
    
    public DeflaterDynHeader(DeflateEncodeTable litLenTable, DeflateEncodeTable distLenTable) {
        this.litLenTable = litLenTable;
        this.distLenTable = distLenTable;
    }
    
    public final void write(final BitOutputStream out) throws IOException {
        
        int hlit = litLenTable.bit_lengths.length;
        while (litLenTable.bit_lengths[--hlit] == 0 && hlit > 256) {}
        out.writeBits(hlit - 256, 5); // 5 bits - the number of Literal/Length codes - 257 (257 - 286)
        
        int hdist = distLenTable.bit_lengths.length;
        while (distLenTable.bit_lengths[--hdist] == 0 && hdist > 0) {}
        out.writeBits(hdist, 5);  // 5 bits - the number of Distance codes - 1         (1 - 32)
        
        hlit++;
        hdist++;
        
        final byte[] len = new byte[hlit + hdist];
        System.arraycopy(litLenTable.bit_lengths, 0, len, 0, hlit);
        System.arraycopy(distLenTable.bit_lengths, 0, len, hlit, hdist);

        new DeflateLengthsTable(len).write(out);
    }
}
