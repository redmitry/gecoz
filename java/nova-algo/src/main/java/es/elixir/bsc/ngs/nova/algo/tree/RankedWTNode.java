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

import es.elixir.bsc.ngs.nova.io.AbstractBitStream;
import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * @author Dmitry Repchevsky
 */

public class RankedWTNode extends AbstractBitStream implements WaveletTreeNode {

    public RankedWTNode(ByteBuffer buf, long size) {
        super(buf.slice(), size);
        super.buf.limit(bytes(size));
        buf.position(buf.position() + super.buf.limit());
    }

    /**
     * Creates the new Wavelet Tree Node.
     * 
     * @param size - the size of the node in bits;
     */
    public RankedWTNode(long size) {
        super(ByteBuffer.allocate(bytes(size)), size);
    }

    /**
     * Calculate the real size of the node in bytes
     * 
     * @param len the length of the bit vector in bits
     * 
     * @return the size of the ranked bit vector structure in bytes
     */
    public final static int bytes(long len) {
        // every 'long' counter replaces the last 'short' counter so 8 - 2 = 6
        final long size = ((len - 1) >>> 16) * 6 + ((len - 1) >>> 9) * 2 + ((len + 7) >>> 3);
        if (size > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("The ranked bit verctor implementation supports ~15G bits");
        }
        return (int)size;
    }
    
    @Override
    public long size() {
        return size;
    }

    /**
     * Reads a bit at the given position.
     * 
     * @param idx
     * @return either '0' or '1'
     */
    @Override
    public byte get(long idx) {
        final long pos = (idx >>> 3) + (idx >>> 9) * 2 + (idx >>> 16) * 6;
        return (byte)((buf.get((int)pos) >>> (idx & 7)) & 0x01);
    }

    @Override
    public void align() throws IOException {
        bits_left &= 0b11111000;
    }

    /**
     * Count the number of 'one' bits up to the position (rank).
     * 
     * @param idx
     * @return 
     */
    @Override
    public long count(long idx) {

        long count = 0;

        final long nlidx = idx >>> 16; // number of large indexes (every 64 kbits)
        final long nsidx = (idx >>> 9) & 0b1111111; // number of small indexes in a block
        final long spos = nsidx * 66; // position of the byte stream
        
        long lpos = 0;
        if (nlidx > 0) {
            lpos = nlidx * 8454; // large index segment start position (nlidx * (8192 + 8 + 127 * 2))
            count = buf.getLong((int)(lpos - 8));
        }
        
        long bpos = lpos + spos;
        if (nsidx > 0) {
            count += buf.getShort((int)(bpos - 2)) & 0xFFFF;
        }
        
        for (int n = (int)(bpos + ((idx >>> 3) & 0b111000)); bpos < n; bpos += 8) {
            count += Long.bitCount(super.getLong(bpos));
        }

        return count + Long.bitCount(super.getLong(bpos) << (63 - (idx & 63)));
    }

    /**
     * Find the position of the n'th zero bit in the node.
     * 
     * @param n
     * @return 
     */
    public long findZero(long n) {
        return findZero(n, n - 1, size - 1);
    }

    /**
     * Find the position of the n'th one bit in the node.
     * 
     * @param n the number of '1' bit to find
     * 
     * @return the position of n'th one bit in the node.
     */
    public long findOne(long n) {
        return findOne(n, n - 1, size - 1);
    }

    protected long findZero(long n, long lo, long hi) {        
        while (lo < hi) {
            final long clo = lo - count(lo) + 1;
            final long chi = hi - count(hi) + 1;
            if (clo >= chi) { // clo > chi would be definitely a bug
                return n == clo && get(lo) == 0 ? lo : -1;
            }
            final long mid = lo + (long)(((hi - lo) * (double)(n - clo)) / (chi - clo));
            final long cmid = mid - count(mid) + 1;
            if (n < cmid) {
                hi = mid - 1;
            } else if (n > cmid) {
                lo = mid + 1;
            } else if (get(mid) == 0) {
                return mid;
            } else {
                hi = mid - 1;
            }
        }
        if (lo == hi && get(lo) == 0 && n == lo - count(lo) + 1) {
            return lo;
        }

        return -1;
    }

    protected long findOne(long n, long lo, long hi) {
        while (lo < hi) {
            final long clo = count(lo);
            final long chi = count(hi);
            if (clo >= chi) {
                return n == clo && get(lo) > 0 ? lo : -1;
            }
            final long mid = lo + (long)(((hi - lo) * (double)(n - clo)) / (chi - clo));
            final long cmid = count(mid);
            if (n < cmid) {
                hi = mid - 1;
            } else if (n > cmid) {
                lo = mid + 1;
            } else if (get(mid) > 0) {
                return mid;
            } else {
                hi = mid - 1;
            }
        }
        if (lo == hi && get(lo) > 0 && n == count(lo)) {
            return lo;
        }
        return -1;
    }
            
    @Override
    protected long getLong() {

        final int pos = buf.position();
        int nlong = pos - ((pos / 8454) * 6);
        nlong -= ((nlong / 66) << 1);

        if ((nlong & 0x00001FFF) == 0 && nlong > 0) {
            buf.position(pos + 8); // skip big index;
        } else if ((nlong & 0b111111) == 0 && nlong > 0) {
            buf.position(pos + 2); // skip small index;
        }
        
        return super.getLong();
    }
    
    @Override
    protected long getLong(long index) {

        long nlong = index - ((index / 8454) * 6); // 8448 + 8 - 2 (when there is a big index no small one needed)
        nlong -= ((nlong / 66) << 1);               // 64 + 2 bytes for small indexes

        if ((nlong & 0x00001FFF) == 0 && nlong > 0) {
            index += 8; // skip big index;
        } else if ((nlong & 0b111111) == 0 && nlong > 0) {
            index += 2; // skip small index;
        }

        return super.getLong(index);
    }
    
    @Override
    protected void putLong(long value) {

        final int pos = buf.position();
        long nlong = pos - ((pos / 8454) * 6);
        nlong -= ((nlong / 66) << 1);

        if ((nlong & 0x00001FFF) == 0 && nlong > 0) {
            buf.putLong(count((nlong << 3) - 1)); // pos * 8 bits
        } else if ((nlong & 0b111111) == 0 && nlong > 0) {            
            int count = (nlong & 0x1FFF) > 64 ? buf.getShort(pos - 66) & 0xFFFF : 0;
            for (int i = pos - 64; i < pos; i +=8) {
                count += Long.bitCount(buf.getLong(i));
            }
            buf.putShort((short) count);
        }

        super.putLong(value);
    }
}
