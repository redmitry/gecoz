/**
 * *****************************************************************************
 * Copyright (C) 2016 Spanish National Bioinformatics Institute (INB) and
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

import es.elixir.bsc.ngs.nova.io.BitBuffer;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Dmitry Repchevsky
 */

public class DeflateTablesTest {
    
    private final static String TEXT = "My uncle, what a worthy man, " +
                                       "Falling ill like that, and dying; " +
                                       "It summons up respect, one can " +
                                       "Admire it, as if he were trying. " +
                                       "Let us all follow his example! " +
                                       "But, God, what tedium to sample " +
                                       "That sitting by the bed all day, " +
                                       "All night, barely a foot away! " +
                                       "And the hypocrisy, demeaning, " +
                                       "Of cosseting one who's half alive; " +
                                       "Puffing the pillows, you contrive " +
                                       "To bring his medicine unsmiling, " +
                                       "Thinking with a mournful sigh, " +
                                       "'Why the devil can't you die?'";

    @Test
    public void test_code_gens() throws IOException {

        byte[] b = TEXT.getBytes();
        long[] c = new long[256];
        for (int i = 0, n = b.length; i < n; i++) {
            c[b[i] & 0xFF]++;
        }

        DeflateEncodeTable het = new DeflateEncodeTable(c);
        
        // check that all symbols with zero frequency have no bit_lengths code assigned
        for (int i = 0; i < 256; i++) {
            Assert.assertTrue((c[i] == 0) == (het.bit_lengths[i] == 0));
        }

        DeflateLookupTable hlt = new DeflateLookupTable(het.bit_lengths);
        
        // check that assigned symbols are properly decoded
        for (int i = 0; i < 256; i++) {
            if (c[i] > 0) {
                Assert.assertEquals(i, hlt.getSymbol(het.table[i]));
            }
        }
    }
    
    @Test
    public void test() throws IOException {
        byte[] b = TEXT.getBytes();
        
        BitBuffer bit_buffer = new BitBuffer(ByteBuffer.allocate(252));
        
        long[] c = new long[256];
        for (int i = 0, n = b.length; i < n; i++) {
            c[b[i] & 0xFF]++;
        }

        DeflateEncodeTable het = new DeflateEncodeTable(c);
        
        for (int i = 0, n = b.length; i < n; i++) {
            het.putSymbol(b[i] & 0xFF, bit_buffer);
        }
        
        bit_buffer.rewind();

        DeflateLookupTable hlt = new DeflateLookupTable(het.bit_lengths);
        
        for (int i = 0, n = b.length; i < n; i++) {
            int symbol = hlt.getSymbol(bit_buffer);
            Assert.assertEquals(symbol, (b[i] & 0xFF));
        }
    }
    
    @Test
    public void encode_test() {
//        ByteArrayOutputStream out = new ByteArrayOutputStream();
//        for (int i = 0; i < 256; i++) {
//            byte[] b = new byte[i * i + 1];
//            Arrays.fill(b, (byte)i);
//            out.write(b, 0, b.length);
//        }
//        
//        byte[] data = out.toByteArray();
//
//        int[] counts = new int[256];
//        for (int i = 0, n = data.length; i < n; i++) {
//            counts[data[i] & 0xFF]++;
//        }
//
//        DeflateEncodeTable encode = new DeflateEncodeTable(counts);
//        for (int i = 0; i < 256; i++) {
//            if (encode.bit_lengths[i] > 0) {
//                System.out.print(String.format("%3s (%2s) ", i ,encode.bit_lengths[i]));
//                System.out.print("                ".substring(encode.bit_lengths[i]));
//                System.out.print(String.format("%" + encode.bit_lengths[i] + "s", Integer.toBinaryString(encode.table[i])).replace(' ', '0'));
//                System.out.println();
//            }
//        }
    }
    
    @Test
    public void stress_test1() throws IOException {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        for (int i = 0; i < 256; i++) {
            byte[] b = new byte[i * i + 1];
            Arrays.fill(b, (byte)i);
            out.write(b, 0, b.length);
        }
        
        byte[] data = out.toByteArray();

        long[] counts = new long[256];
        for (int i = 0, n = data.length; i < n; i++) {
            counts[data[i] & 0xFF]++;
        }

        DeflateEncodeTable encode = new DeflateEncodeTable(counts);
        DeflateLookupTable decode = new DeflateLookupTable(encode.bit_lengths);
        
        for (int i = 0; i < 256; i++) {
            if (counts[i] > 0) {
                if (i != decode.getSymbol(encode.table[i])) {
                    decode.getSymbol(encode.table[i]);
                }
                Assert.assertEquals(i, decode.getSymbol(encode.table[i]));
            }
        }
    }
    
    @Test
    public void stress_test2() throws IOException {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        for (int i = 0; i < 256; i++) {
            byte[] b = new byte[i * i + 1];
            Arrays.fill(b, (byte)i);
            out.write(b, 0, b.length);
        }
        
        byte[] data = out.toByteArray();
        
        BitBuffer bit_buffer = new BitBuffer(ByteBuffer.allocate(data.length));
        
        long[] counts = new long[256];
        for (int i = 0, n = data.length; i < n; i++) {
            counts[data[i] & 0xFF]++;
        }

        DeflateEncodeTable encode = new DeflateEncodeTable(counts);
        
        for (int i = 0, n = data.length; i < n; i++) {
            encode.putSymbol(data[i] & 0xFF, bit_buffer);
        }
        
        bit_buffer.rewind();

        DeflateLookupTable decode = new DeflateLookupTable(encode.bit_lengths);
        
        for (int i = 0, n = data.length; i < n; i++) {
            int symbol = decode.getSymbol(bit_buffer);
            Assert.assertEquals(data[i] & 0xFF, symbol);
        }
    }
}
