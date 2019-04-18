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

import es.elixir.bsc.ngs.nova.io.BitInputStream;
import java.io.IOException;

/**
 * Deflater Block Header class used for the decompression.
 * 
 * @author Dmitry Repchevsky
 */

public class InflaterBlockHeader extends BlockHeader {

    private static DeflateLookupTable defaultLitLenTree;
    private static DeflateLookupTable defaultDistTree;

    public final DeflateLookupTable litLenTree;
    public final DeflateLookupTable distTree;
    
    public InflaterBlockHeader(final BitInputStream in) throws IOException {
        
        super(in);
        
        switch(btype) {
            case NO_COMPRESSION: litLenTree = null;
                                 distTree = null;

                                 in.skipBits(5);
                                 len = (int)(in.readBits(16) & 0xFFFF);
                                 final int nlen = (int)(in.readBits(16) & 0xFFFF);
                                 if (nlen != (len ^ 0xFFFF)) {
                                    throw new IOException("broken uncompressed block");
                                 }
                                 return;
            case HUFF_FIXED:     litLenTree = getDefaultLitLenTree();
                                 distTree = getDefaultDistTree();
                                 break;
            case HUFF_DYNAMIC:   InflaterDynHeader header = new InflaterDynHeader(in);
                                 litLenTree = header.litLenTree;
                                 distTree = header.distTree;
                                 break;
            default: throw new IOException("invalid btype");
        }
    }

    public static final synchronized DeflateLookupTable getDefaultLitLenTree() {
        if (defaultLitLenTree == null) {
            defaultLitLenTree = new DeflateLookupTable(getDefaultLitLen());
        }
        
        return defaultLitLenTree;
    }
    
    public static final synchronized DeflateLookupTable getDefaultDistTree() {
        if (defaultDistTree == null) {
            defaultDistTree = new DeflateLookupTable(getDefaultDist());
        }
        
        return defaultDistTree;
    }
}
