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

import java.nio.charset.StandardCharsets;

/**
 * @author Dmitry Repchevsky
 */

public class Sequence {
    final static byte[] SEQ = {'=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
    
    private final int l_seq;
    private final byte[] seq; // compact sequence form (2 letters per byte)
    
    public Sequence(final String sequence) {
        l_seq = sequence.length();
        seq = Sequence.parse(sequence);
    }
    
    protected Sequence(final byte[] sequence, final int l_seq) {
        seq = sequence;
        this.l_seq = l_seq;
    }

    @Override
    public String toString() {
        return seq != null ? parse(seq, l_seq) : "";
    }

    public static byte[] parseToBytes(
            final byte[] sequence, final int l_seq) {

        final byte[] arr = new byte[l_seq];
        for (int i = 0, idx = 0, n = l_seq >> 1; i < n; i++) {
            arr[idx++] = SEQ[sequence[i] >> 4 & 0x0F];
            arr[idx++] = SEQ[sequence[i] & 0x0F];
        }
        if ((l_seq & 1) != 0) {
           arr[l_seq - 1] = SEQ[sequence[sequence.length - 1] >> 4 & 0x0F];
        }
        return arr;
    }

    /**
     * <p>
     * Decode binary encoded sequence (SAM v1.6 4.2.3).
     * </p>
     * 
     * @param sequence four bit encoded sequence
     * @param l_seq sequence length
     * 
     * @return decoded SAM sequence (SEQ)
     */
    public static String parse(final byte[] sequence, final int l_seq) {
        return new String(parseToBytes(sequence, l_seq), StandardCharsets.ISO_8859_1);
    }
    
    /**
     * <p>
     * Encode SAM sequence into the binary form (SAM v1.6 4.2.3).
     * </p>
     * 
     * @param sequence SAM sequence to encode
     * 
     * @return encoded in 4-bit values sequence
     */
    public static byte[] parse(final String sequence) {

        final int l_seq = sequence.length();
        final byte[] seq = new byte[(l_seq + 1) >> 1];
        
        final int n = l_seq & 0xFFFFFFFE; // make n even
        for (int i = 0, j = 0; j < n; i++) {
            final byte a = Sequence.getByte(sequence.charAt(j++));
            final byte b = Sequence.getByte(sequence.charAt(j++));
            seq[i] = (byte)(a << 4 | b);
        }
        
        if (l_seq != n) {
            // set the last (odd) character
            seq[seq.length - 1] = (byte)(Sequence.getByte(sequence.charAt(n)) << 4);
        }
        
        return seq;
    }
    
    private static byte getByte(final char ch) {
        byte i;
        switch(Character.toUpperCase(ch)) {
            case '=': i = 0x00; break;
            case 'A': i = 0x01; break;
            case 'C': i = 0x02; break;
            case 'M': i = 0x03; break;
            case 'G': i = 0x04; break;
            case 'R': i = 0x05; break;
            case 'S': i = 0x06; break;
            case 'V': i = 0x07; break;
            case 'T': i = 0x08; break;
            case 'W': i = 0x09; break;
            case 'Y': i = 0x0A; break;
            case 'H': i = 0x0B; break;
            case 'K': i = 0x0C; break;
            case 'D': i = 0x0D; break;
            case 'B': i = 0x0E; break;
            case 'N': i = 0x0F; break;
            default: i = Byte.MIN_VALUE;
        }
        return i;
    }
}
