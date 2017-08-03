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
import java.util.Arrays;

/**
 * @author Dmitry Repchevsky
 */

public class BlockHeader {

    private static DeflateLookupTable defaultLitLenTree;
    private static DeflateLookupTable defaultDistTree;

    public final static int NO_COMPRESSION = 0x00;
    public final static int HUFF_FIXED = 0x01;
    public final static int HUFF_DYNAMIC = 0x02;

    public final boolean bfinal;
    public final int btype;
    public int len;

    public final DeflateLookupTable litLenTree;
    public final DeflateLookupTable distTree;
    
    public BlockHeader(BitInputStream in) throws IOException {
        final int block_header = (int) ((in.readBits(3)) & 0b111);
        
        bfinal = (block_header & 0x01) != 0;
        btype = block_header >>> 1;
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
            case HUFF_FIXED:     litLenTree = BlockHeader.getDefaultLitLenTree();
                                 distTree = BlockHeader.getDefaultDistTree();
                                 break;
            case HUFF_DYNAMIC:   InflaterDynHeader header = new InflaterDynHeader(in);
                                 litLenTree = header.litLenTree;
                                 distTree = header.distTree;
                                 break;
            default: throw new IOException("invalid btype");
        }
    }
    
    public static synchronized DeflateLookupTable getDefaultLitLenTree() {
        if (defaultLitLenTree == null) {
            final byte[] cl = new byte[288];
            Arrays.fill(cl, 0, 144, (byte)8);
            Arrays.fill(cl, 144, 256, (byte)9);
            Arrays.fill(cl, 256, 280, (byte)7);
            Arrays.fill(cl, 280, 288, (byte)8);
            defaultLitLenTree = new DeflateLookupTable(cl);
        }
        
        return defaultLitLenTree;
    }
    
    public static synchronized DeflateLookupTable getDefaultDistTree() {
        if (defaultDistTree == null) {
            final byte[] cl = new byte[32];
            Arrays.fill(cl, 0, 32, (byte)5);
            defaultDistTree = new DeflateLookupTable(cl);
        }
        
        return defaultDistTree;
    }
}
