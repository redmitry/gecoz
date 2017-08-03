package es.elixir.bsc.ngs.nova.algo.deflate;

import es.elixir.bsc.ngs.nova.io.BitInputStream;
import es.elixir.bsc.ngs.nova.io.BitOutputStream;
import java.io.IOException;

/**
 * @author Dmitry Repchevsky
 */

public class DeflateLengthsTable {
    
    private static final byte[] CL_ORDER = 
        { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };

    public final byte[] d_tree;
    
    public DeflateLengthsTable(BitInputStream in, int len) throws IOException {
        this.d_tree = new byte[len];
        
        final int hclen = (int) ((in.readBits(4) & 0b1111) + 4); // 4 bits - the number of Code Length codes - 4     (4 - 19)
        
        final byte[] l_tree = new byte[19];
        for (int i = 0; i < hclen; i++) {
            l_tree[CL_ORDER[i]] = (byte)(in.readBits(3) & 0b111);
        }
        
        DeflateLookupTable table = new DeflateLookupTable(l_tree);
        byte symbol = 0;

        for (int i = 0, n = d_tree.length; i < n;) {
            final byte code = (byte)table.getSymbol(in);
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
    }
            
    public DeflateLengthsTable(byte[] d_tree) {
        this.d_tree = d_tree;
    }
    
    public final void write(final BitOutputStream out) throws IOException {
        
        final long[] counts = new long[19];
        final int hclen = counts(d_tree, counts);
        
        out.writeBits(hclen - 3, 4);

        DeflateEncodeTable table = new DeflateEncodeTable(counts);
        
        for (int i = 0; i <= hclen; i++) {
            out.writeBits(table.bit_lengths[CL_ORDER[i]], 3);
        }
        
        for (int i = 0, len = 0, count = 0, n = d_tree.length - 1; i <= n; i++) {
            if (len != d_tree[i] || i == n) {
                while (count >= 3) {
                    if (len != 0) {
                        table.putSymbol(16, out);
                        count -= 3;
                        out.writeBits(Math.min(count, 3), 2);
                        count -= 3;
                    } else if (count <= 10) {
                        table.putSymbol(17, out);
                        count -= 3;
                        out.writeBits(Math.min(count, 7), 3);
                        count -= 7;
                    } else {
                        table.putSymbol(18, out);
                        count -= 11;
                        out.writeBits(Math.min(count, 127), 7);
                        count -= 127;
                    }
                }
                while (count-- > 0) {
                    table.putSymbol(len, out);
                }
                len = d_tree[i];
                table.putSymbol(len, out);
                count = 0;
            } else {
                count++;
            }
        }
    }
    
    /**
     * <p>
     * Calculates the length of the table.
     * </p>
     * 
     * @param bit_lengths - Deflate bit lengths table.
     * 
     * @return the length (in bits) of the table.
     */
    public static int length(byte[] bit_lengths) {
        
        final long[] counts = new long[19];
        final int hclen = counts(bit_lengths, counts);
        
        int bits = 7 + hclen * 3;
        
        DeflateEncodeTable table = new DeflateEncodeTable(counts);
        
        for (int i = 0, len = 0, count = 0, n = bit_lengths.length - 1; i <= n; i++) {
            if (len != bit_lengths[i] || i == n) {
                while (count >= 3) {
                    if (len != 0) {
                        bits += table.bit_lengths[16] + 2;
                        count -= 6;
                    } else if (count <= 10) {
                        bits += table.bit_lengths[17] + 3;
                        count -= 10;
                    } else {
                        bits += table.bit_lengths[18] + 7;
                        count -= 138;
                    }
                }
                while (count-- > 0) {
                    bits += table.bit_lengths[len];
                }
                len = bit_lengths[i];
                bits += table.bit_lengths[len];
                count = 0;
            } else {
                count++;
            }
        }
        
        return bits;
    }

    private static int counts(byte[] bit_lengths, long[] counts) {
        for (int i = 0, len = 0, count = 0, n = bit_lengths.length - 1; i <= n; i++) {
            if (len != bit_lengths[i] || i == n) {
                while (count >= 3) {
                    if (len != 0) {
                        counts[16]++;
                        count -= 6;
                    } else if (count <= 10) {
                        counts[17]++;
                        count -= 10;
                    } else {
                        counts[18]++;
                        count -= 138;
                    }
                }
                while (count-- > 0) {
                    counts[len]++;
                }
                len = bit_lengths[i];
                counts[len]++;
                count = 0;
            } else {
                count++;
            }
        }
        
        int hclen = 18;
        do {
            if (counts[CL_ORDER[hclen]] > 0) {
                break;
            }            
        } while (--hclen >= 0);
        
        return hclen;
    }

}
