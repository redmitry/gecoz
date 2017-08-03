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

import java.io.IOException;
import java.util.Arrays;
import java.util.zip.Checksum;

/**
 * @author Dmitry Repchevsky
 */
class InflaterOutput {
    
    private final static int MASK = 0x7FFF;
    
    private final Checksum hash;
    private final byte[] buffer;
    private int start, end;
    private int count;
    
    InflaterOutput(Checksum hash) {
        this.hash = hash;
        buffer = new byte[32768];
    }
    
    public void reset() {
        count = 0;
        start = 0;
        end = 0;
        if (hash != null) {
            hash.reset();
        }
    }

    /**
     * @return the number of bytes already read
     */
    public int count() {
        return count;
    }

    /**
     * @return the number of bytes available to read
     */
    public int available() {
        return (end - start) & MASK;
    }

    public void push(int i) {
        buffer[end++] = (byte)i;
        end &= MASK;
    }
    
    public void inflate(int dist, int len) {
        if (dist < len) {
            overlap(dist, len);
        } else if (dist <= end) {
            final int pos = end - dist;
            final int size = buffer.length - end;
            if (size >= len) {
                System.arraycopy(buffer, pos, buffer, end, len);
                end = (end + len) & MASK;
            } else {
                System.arraycopy(buffer, pos, buffer, end, size);
                end = len - size;
                System.arraycopy(buffer, buffer.length - dist, buffer, 0, end);
            }
        } else {
            final int size = dist - end;
            if (size > len) {
                System.arraycopy(buffer, buffer.length - size, buffer, end, len);
                end = (end + len) & MASK;
            } else {
                System.arraycopy(buffer, buffer.length - size, buffer, end, size);
                end +=size;
                inflate(end, len - size);
            }
        }
    }
            
    private void slow(int pos, int len) {
        while(len-- > 0) {
            buffer[end++] = buffer[pos++];
            end &= MASK;
            pos &= MASK;
        }
    }
    
    private void overlap(int dist, int len) {
        int pos = (end - dist) & MASK;
        
        if (dist == 1) {
            // repeat the last symbol 'len' times
            final byte ch = buffer[pos];
            final int size = buffer.length - end;
            if (size > len) {
                Arrays.fill(buffer, end, end + len, ch);
                end = (end + len) & MASK;
            } else {
                Arrays.fill(buffer, end, buffer.length, ch);
                end = len - size;
                Arrays.fill(buffer, 0, end, ch);
            }
        } else {
            slow(pos, len);
        }
    }
    
    public int read() {
        if (end == start) {
            return -1;
        }
        
        final int b = buffer[start] & 0xFF;
        if (hash != null) {
            hash.update(b);
        }
        count++;
        start++;
        start &= MASK;
        return b;
    }

    public final long skip(long nbytes) throws IOException {
        nbytes = Math.min(nbytes, available());
        if (hash != null) {
            final int size = buffer.length - start;
            if (size >= nbytes) {
                hash.update(buffer, start, (int)nbytes);
            } else {
                hash.update(buffer, start, (int)size);
                hash.update(buffer, 0, (int)(nbytes - size));
                
            }
        }
        count += nbytes;
        start += nbytes;
        start &= MASK;
        
        return nbytes;
    }
    
    public int read(byte[] buf) {
        return read(buf, 0, buf.length);
    }
    
    public int read(byte[] buf, int off, int len) {
        if (end == start) {
            return 0;
        } else if (end > start) {
            final int size = Math.min(len, end - start);
            System.arraycopy(buffer, start, buf, off, size);
            if (hash != null) {
                hash.update(buf, off, size);
            }
            start = (start + size) & MASK;
            count += size;
            return size;
        } else {
            int size = buffer.length - start;
            if (size >= len) {
                System.arraycopy(buffer, start, buf, off, len);
                if (hash != null) {
                    hash.update(buf, off, len);
                }
                start = (start + len) & MASK;
                count += len;
                return len;
            } else {
                System.arraycopy(buffer, start, buf, off, size);
                if (hash != null) {
                    hash.update(buf, off, size);
                }
                start = Math.min(end, len - size);
                System.arraycopy(buffer, 0, buf, off + size, start);
                if (hash != null) {
                    hash.update(buf, off + size, start);
                }
                size += start;
                count += size;
                return size;
            }   
        }
    }
}
