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

package es.elixir.bsc.ngs.nova.algo.tree;

import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * The Huffman Shaped Wavelet Tree implementation.
 * The implementation uses a table of nodes indexed via Huffman codes.
 * 
 * @author Dmitry Repchevsky
 */

public class HuffmanShapedWaveletTree {

    public final HSWTShape shape;
    private final RankedWTNode[] nodes;

    /**
     * Creates Huffman Shaped Wavelet Tree from a byte array.
     * 
     * @param s
     * @throws IOException 
     */
    public HuffmanShapedWaveletTree(byte[] s) throws IOException {
        this(ByteBuffer.wrap(s));
    }

    /**
     * Creates Huffman Shaped Wavelet Tree from a byte buffer.
     * 
     * @param src
     * @throws IOException 
     */
    public HuffmanShapedWaveletTree(ByteBuffer src) throws IOException {
        this(new ByteBufferDataSource(src));
    }

    /**
     * Creates Huffman Shaped Wavelet Tree from a byte buffer.
     * 
     * @param src
     * @throws IOException 
     */
    public HuffmanShapedWaveletTree(DataSource src) throws IOException {
        this(src, null);
    }
    
    /**
     * Creates Huffman Shaped Wavelet Tree into provided byte buffer.
     * 
     * @param src buffer that contains string data
     * @param dst buffer to be used to keep WT vectors
     * 
     * @throws IOException 
     */
    public HuffmanShapedWaveletTree(DataSource src, ByteBuffer dst) throws IOException {
        this(new HSWTShape(count(src)), src, dst);
    }
    
    /**
     * Creates Huffman Shaped Wavelet Tree with predefined shape. 
     * Note that the shape must be created from the provided source.
     * 
     * @param shape the shape of the tree
     * @param src
     * @param dst
     * @throws IOException 
     */

    private HuffmanShapedWaveletTree(HSWTShape shape, DataSource src, ByteBuffer dst) throws IOException {
        this.shape = shape;

        nodes = new RankedWTNode[256];

        int[] lengths = new int[256]; // the size of bit verctors (in bits)

        for (int i = 0; i < 256; i++) {
            if (shape.counts[i] > 0) {
                final int code = shape.encode.table[i];
                for (int j = 0, n = shape.encode.bit_lengths[i]; j < n; j++) {
                    int idx = code & (0x0000FFFF >>> (16 - j));
                    idx |= (0x8000 >>> (15 - j));
                    idx = shape.decode.getSymbol(idx);
                    lengths[idx] += shape.counts[i];
                }
            }
        }

        if (dst != null) {
            mapNodes(dst, lengths, 1);
        } else {
            for (int i = 0; i < 256; i++) {
                if (lengths[i] > 0) {
                    nodes[i] = new RankedWTNode(lengths[i]);
                }
            }            
        }
        
        fill(src);
    }

    private void fill(DataSource in) throws IOException {
        for (int i = 0, n = in.length(); i < n; i++) {
            final int symbol = in.get(i) & 0xFF;
            final int code = shape.encode.table[symbol];

            for (int j = 0, m = shape.encode.bit_lengths[symbol]; j < m; j++) {
                int idx = code & (0x0000FFFF >>> (16 - j));
                idx |= (0x8000 >>> (15 - j));
                idx = shape.decode.getSymbol(idx);

                nodes[idx].put((code >>> j) & 0x01);
            }            
        }

        for (int i = 0; i < 256; i++) {
            if (nodes[i] != null) {
                nodes[i].flush();
            }
        }        
    }

    /**
     * Creates Huffman Shaped Wavelet Tree with predefined shape and its
     * serialized form.
     * 
     * @param shape the shape of the tree
     * @param in byte buffer that contains a serialized HSWT.
     * 
     * @throws IOException 
     */
    private HuffmanShapedWaveletTree(HSWTShape shape, ByteBuffer in) throws IOException {
        this.shape = shape;
        
        nodes = new RankedWTNode[256];
                
        mapNodes(in, shape.length, 1);
    }

    private void mapNodes(ByteBuffer in, int[] lengths, int code) {

        int idx = shape.decode.getSymbol(code);

        final int level = Integer.numberOfLeadingZeros(code) - 1;
        if (31 - level > shape.encode.bit_lengths[idx]) {
            return; // the leaf
        }

        code |= Integer.MIN_VALUE >>> level;
        
        if (nodes[idx] == null) {
            nodes[idx] = new RankedWTNode(in, lengths[idx]);
            
            mapNodes(in, lengths, code & (0xBFFFFFFF >> level));
            mapNodes(in, lengths, code | (0x40000000 >> level));
        }
    }
    
    /**
     * <p>
     * Reconstructs the tree from a plain ByteBuffer where all nodes go one 
     * after another (byte aligned) in a Huffman codes order.
     * </p>
     * <p>
     * All the nodes data is a view of the provided ByteBuffer.
     * </p>
     * 
     * @param in - the ByteBuffer that contains all the tree nodes
     * @param length - the size of the node (the root node equals text size)
     * @param code - the node's Huffman code
     */
    private void mapNodes(ByteBuffer in, long length, int code) {

        int idx = shape.decode.getSymbol(code);
        final int level = Integer.numberOfLeadingZeros(code) - 1;
        
        if (31 - level > shape.encode.bit_lengths[idx]) {
            return; // the leaf
        }
        
        code |= Integer.MIN_VALUE >>> level;
        
        if (nodes[idx] == null) {
            nodes[idx] = new RankedWTNode(in, length);
            
            final long bits = nodes[idx].count(length - 1);        

            mapNodes(in, length - bits, code & (0xBFFFFFFF >> level)); // left
            mapNodes(in, bits, code | (0x40000000 >> level));          // right
        }
    }
    
    public void write(ByteBuffer out) throws IOException {
        writeNodes(out, 0, 0);
    }

    private void writeNodes(ByteBuffer out, int code, int level) throws IOException {

        int symbol = shape.decode.getSymbol(code, level);
        if (symbol >= 0) {
            return;
        }

        int idx = code | (0x8000 >>> (15 - level));
        idx = shape.decode.getSymbol(idx);
        
        nodes[idx].write(out);
        
        writeNodes(out, code, level + 1); // left
        writeNodes(out, code | (1 << level), level + 1); // right
    }

    
    /**
     * Counts symbols up to the position.
     * 
     * @param symbol
     * @param pos
     * 
     * @return the number of occurrences of the symbol up to the specified position.
     */
    public long occ(int symbol, long pos) {
        if (shape.encode.bit_lengths[symbol] == 0) {
            return -1;
        }

        final int code = shape.encode.table[symbol];
        for (int i = 0, n = shape.encode.bit_lengths[symbol]; i < n && pos >= 0; i++) {
            int idx = (code & (0x0000FFFF >>> (16 - i)));
            idx |= (0x8000 >>> (15 - i));
            idx = shape.decode.getSymbol(idx);

            long bits = nodes[idx].count(pos);
            if (((code >>> i) & 0x01) == 0) {
                pos -= bits;
            } else {
                pos = bits - 1;
            }

        }
        return pos;
    }
    
    public long getRank(long pos) {
        int idx = shape.decode.getSymbol(1);
        
        for (int i = 0, code = 0; i < shape.encode.bit_lengths[idx]; i++) {
            final int bit = nodes[idx].get(pos);
            final long bits = nodes[idx].count(pos);
            pos = bit == 0 ? pos - bits : bits - 1;

            code |= bit << i;
            idx = shape.decode.getSymbol(code | (0x8000 >>> (14 - i)));
        }
        
        return pos;      
    }
    
    public int getSymbol(long pos) {

        int idx = shape.decode.getSymbol(1);
        
        for (int i = 0, code = 0; i < shape.encode.bit_lengths[idx]; i++) {
            final int bit = nodes[idx].get(pos);
            final long bits = nodes[idx].count(pos);
            pos = bit == 0 ? pos - bits : bits - 1;

            code |= bit << i;
            idx = shape.decode.getSymbol(code | (0x8000 >>> (14 - i)));
        }
        
        return idx;
    }
    
    public long getRS(long pos) {
        
        int idx = shape.decode.getSymbol(1);
        
        for (int i = 0, code = 0; i < shape.encode.bit_lengths[idx]; i++) {
            final int bit = nodes[idx].get(pos);
            final long bits = nodes[idx].count(pos);
            pos = bit == 0 ? pos - bits : bits - 1;

            code |= bit << i;
            idx = shape.decode.getSymbol(code | (0x8000 >>> (14 - i)));
        }
        
        return ((long)pos << 32) | idx;
    }

    /**
     * Creates the HSWT from the serialized (raw) form.
     * 
     * @param shape the Huffman's "shape" of the tree
     * @param in byte buffer where the HSWT data (nodes) are stored
     * 
     * @return created HSWT
     * 
     * @throws IOException 
     */
    public static HuffmanShapedWaveletTree read(HSWTShape shape, ByteBuffer in) throws IOException {
        return new HuffmanShapedWaveletTree(shape, in);
    }
    
    public static HuffmanShapedWaveletTree write(HSWTShape shape, DataSource src, ByteBuffer dst) throws IOException {
        return new HuffmanShapedWaveletTree(shape, src, dst);
    }
    
    public static long[] count(DataSource src) {
        long[] counts = new long[256];
        for (int i = 0, n = src.length(); i < n; i++) {
            counts[src.get(i) & 0xFF]++;
        }
        return counts;
    }
    
    public static interface DataSource {
        byte get(int idx);
        int length();
    }
    
    public final static class ByteBufferDataSource implements DataSource {
        
        private final ByteBuffer buf;
        
        public ByteBufferDataSource(final ByteBuffer buf) {
            this.buf = buf;
        }
        
        @Override
        public final byte get(int idx) {
            return buf.get(idx);
        }
        
        @Override
        public final int length() {
            return buf.limit();
        }
    }
}
