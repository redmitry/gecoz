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
import java.util.zip.Checksum;

/**
 * @author Dmitry Repchevsky
 */

public class Inflater {

    // length codes (RFC 1951 3.2.5)
    private static final int BLC_LENS[] =  { 
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258
    };
  
    // extra bits for the lengths
    private static final int BLC_LENS_EXT_BITS[] = { 
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 
        2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0
    };
  
    // block distances ( 0 - 29 )
    private static final int BLC_DIST[] = {
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
        8193, 12289, 16385, 24577
    };
  
    // extra bits for distances offsets
    private static final int BLC_DIST_EXT_BITS[] = {
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13
    };
  
    private final BitInputStream in;
    private final InflaterOutput out;
    private BlockHeader header;
    private boolean end_of_stream;
    
    public Inflater(BitInputStream in, Checksum hash) throws IOException {
        this.in = in;
        out = new InflaterOutput(hash);
    }

    public void reset() {
        header = null;
        end_of_stream = false;
        out.reset();
    }

    /**
     * Checks whether the inflator stream is finished.
     * 
     * @return true if the last inflator block has been completely read.
     * @throws IOException 
     */
    public boolean finished() throws IOException {
        if (out.available() > 0) {
            return false;
        }
        
        if (end_of_stream) {
            return true;
        }
        
        if (header == null) {
            header = new BlockHeader(in);
        }
            
        return inflate(0, 1) == 0;

    }

    /**
     * @return the number of bytes read from the inflater
     */
    public int count() {
        return out.count();
    }

    public final long skip(final long n) throws IOException {
        int available = out.available();
        if (available == 0) {
            if (end_of_stream) {
                return -1;
            } else if (n == 0) {
                return 0;
            } else {
                if (header == null) {
                    header = new BlockHeader(in);
                }

                if (header.btype == BlockHeader.NO_COMPRESSION) {
                    if (header.len > n) {
                        header.len -= n;
                        return n;
                    } else {
                        if (header.bfinal) {
                            end_of_stream = true;
                        } else {
                            header = null;
                        }
                        return header.len;
                    }
                }

                available = inflate(0, (int)Math.min(n, 32510));
                if (available == 0) {
                    return skip(n);
                }
            }
        }

        return out.skip(n);
    }
    
    public int read() throws IOException {
        int available = out.available();
        if (available == 0) {
            if (end_of_stream) {
                return -1;
            }

            if (header == null) {
                header = new BlockHeader(in);
            }
            
            if (header.btype == BlockHeader.NO_COMPRESSION) {
                if (--header.len == 0) {
                    if (header.bfinal) {
                        end_of_stream = true;
                    } else {
                        header = null;
                    }
                    return read();
                }
            } else if (inflate(0, 16) == 0) {
                return read();
            }
        }
        
        return out.read();
    }
    
    public int read(byte[] buf) throws IOException {
        return read(buf, 0, buf.length);
    }
    
    public int read(byte[] buf, int off, int len) throws IOException {
        int available = out.available();
        if (available < len) {
            if (end_of_stream) {
                if (available == 0) {
                    return -1;
                }
                return out.read(buf, off, available);
            }

            if (header == null) {
                header = new BlockHeader(in);
            }
        
            if (header.btype == BlockHeader.NO_COMPRESSION) {
                len = len - available;
                
                // read the rests from previous block;
                while(true) {
                    final int read = out.read(buf, off, available);
                    if (read <= 0) {
                        break;
                    }
                    off += read;
                    available -= read;   
                }
                if (header.len > len) {
                    header.len -= len;
                } else {
                    len = header.len;
                    if (header.bfinal) {
                        end_of_stream = true;
                    } else {
                        header = null;
                    }
                }
                
                //return in.read(buf, off, len);
                
                in.align();
                for (int i = 0; i < len; i++) {
                    buf[off++] = (byte)(in.readBits(8) & 0xFF);
                }
                return len;
                
            }
            // ensure we do not exceed the buf.size (32768 - 258 = 32510)
            return out.read(buf, off, inflate(available, Math.min(len, 32510)));
        }
        return out.read(buf, off, len);
    }

    /**
     * Fills the buffer with decompressed data
     *
     * @param available the unread data that is already available in the buffer
     * @param len the data we would like to have in the buffer
     * 
     * @return the length of data we can safely read from the buffer 
     * ( available &lt;= length &lt;= len )
     * 
     * @throws IOException 
     */
    private int inflate(int available, int len) throws IOException {
        while (available <= len) {
            final int length = inflate();
            if (length == 0) {
                if (header.bfinal) {
                    end_of_stream = true;
                } else {
                    header = null;
                }
                break;
            }
            available += length;
        }
        return Math.min(available, len);
    }
    
    private int inflate() throws IOException {
        final int symbol = header.litLenTree.getSymbol(in);
        if (symbol <= 255) {
            out.push(symbol);
            return 1;
        } else if (symbol == 256) {
            return 0;
        }
        final int ext_len_bits = BLC_LENS_EXT_BITS[symbol - 257];
        final int length = (int) (ext_len_bits == 0 ? BLC_LENS[symbol - 257] : 
                BLC_LENS[symbol - 257] + (in.readBits(ext_len_bits) & (0xFFFFFFFF >>> (32 - ext_len_bits))));

        final int dist = header.distTree.getSymbol(in);
        final int ext_dist_bits = BLC_DIST_EXT_BITS[dist];
        final int distance = (int) (ext_dist_bits == 0 ? BLC_DIST[dist] :
                BLC_DIST[dist] + (in.readBits(ext_dist_bits) & (0xFFFFFFFF >>> (32 - ext_dist_bits))));

        out.inflate(distance, length);
        return length;        
    }
}
