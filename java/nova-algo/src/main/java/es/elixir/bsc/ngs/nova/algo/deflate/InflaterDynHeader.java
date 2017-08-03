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
    
    private static final byte[] CL_ORDER = 
        { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };

    public InflaterDynHeader(final BitInputStream in) throws IOException {
        final int hlit = (int) ((in.readBits(5) & 0b11111) + 257); // 5 bits - the number of Literal/Length codes - 257 (257 - 286)
        final int hdist = (int) ((in.readBits(5) & 0b11111) + 1); // 5 bits - the number of Distance codes - 1        (1 - 32)
        final int hclen = (int) ((in.readBits(4) & 0b1111) + 4); // 4 bits - the number of Code Length codes - 4     (4 - 19)
        
        final byte[] l_tree = new byte[19];
        for (int i = 0; i < hclen; i++) {
            l_tree[CL_ORDER[i]] = (byte)(in.readBits(3) & 0b111);
        }
        
        DeflateLookupTable literalTree = new DeflateLookupTable(l_tree);
        byte symbol = 0;
        final byte[] d_tree = new byte[hlit + hdist];
        for (int i = 0, n = d_tree.length; i < n;) {
            final byte code = (byte)literalTree.getSymbol(in);
            if (code <= 15) {
                d_tree[i++] = symbol = code;
            } else if (code == 16) {
                // repeate symbol 3-6 times
                final int rep = (int) ((in.readBits(2) & 0b11) + 3);
                for (int j = 0; j < rep; j++, i++) {
                    d_tree[i] = symbol;
                }
            } else if (code == 17) {
                // repeate zeros 3 - 10 times (just skip the index)
                final int rep = (int) ((in.readBits(3) & 0b111) + 3);
                i+= rep;
            } else if (code == 18) {
                // repeate zeros 11 - 138 times (just skip the index)
                final int rep = (int) ((in.readBits(7) & 0b1111111) + 11);
                i+= rep;
            }
        }
        
        final byte[] litLenCodes = new byte[hlit];
        System.arraycopy(d_tree, 0, litLenCodes, 0, hlit);
        litLenTree = new DeflateLookupTable(litLenCodes);
        
        final byte[] distLens = new byte[hdist];
        System.arraycopy(d_tree, hlit, distLens, 0, hdist);
        distTree = new DeflateLookupTable(distLens);
    }
}
