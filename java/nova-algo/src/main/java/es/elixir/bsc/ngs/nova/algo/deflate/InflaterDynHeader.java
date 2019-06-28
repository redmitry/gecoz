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

package es.elixir.bsc.ngs.nova.algo.deflate;

import es.elixir.bsc.ngs.nova.io.BitInputStream;
import java.io.IOException;

/**
 * @author Dmitry Repchevsky
 */

public class InflaterDynHeader {

    public final DeflateLookupTable litLenTree;
    public final DeflateLookupTable distTree;

    public InflaterDynHeader(final BitInputStream in) throws IOException {
        final int hlit = (int) ((in.readBits(5) & 0b11111) + 257); // 5 bits - the number of Literal/Length codes - 257 (257 - 286)
        final int hdist = (int) ((in.readBits(5) & 0b11111) + 1);  // 5 bits - the number of Distance codes - 1         (1 - 32)
        
        final DeflateLengthsTable lengthsTable = new DeflateLengthsTable(in, hlit + hdist);

        final byte[] litLenCodes = new byte[hlit];
        System.arraycopy(lengthsTable.d_tree, 0, litLenCodes, 0, hlit);
        litLenTree = new DeflateLookupTable(litLenCodes);
        
        final byte[] distLenCodes = new byte[hdist];
        System.arraycopy(lengthsTable.d_tree, hlit, distLenCodes, 0, hdist);
        distTree = new DeflateLookupTable(distLenCodes);
    }
}
