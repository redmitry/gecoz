/**
 * *****************************************************************************
 * Copyright (C) 2017 Spanish National Bioinformatics Institute (INB) and
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

package es.elixir.bsc.ngs.nova.gecoz;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.zip.DataFormatException;

/**
 * @author Dmitry Repchevsky
 */

public class GecozRefBlockHeader {
    public final static String MAGIC = "GecozBWT";
    
    public final byte version = 1;
    public final long size; // the block size
    public final long len;  // the length of the generalized string
    public final String[] headers;
    
    /**
     * Block header for the BWT which was generated from one or more sequences.
     * 
     * @param headers a map that has sequences´ headers and their position in the chain
     * @param size the size of the BWT block
     * @param len
     */
    public GecozRefBlockHeader(String[] headers, long size, long len) {
        this.headers = headers;

        this.size = size;
        this.len = len;
    }
    
    public GecozRefBlockHeader(InputStream in) throws IOException, DataFormatException {
        DataInputStream dis = new DataInputStream(in);

        if (dis.readLong() != 0x4765636F7A425754L |
            dis.readByte() != version) { // "GecozBWT" v1
        }
        
        size = Long.reverseBytes(dis.readLong());
        len = Long.reverseBytes(dis.readLong());
        
        ArrayList<String> list = new ArrayList<>();
        // the list of headers separated by 0x00 and terminated by double zero
        int ch;
        while((ch = dis.read()) > 0) {
            StringBuilder sb = new StringBuilder();
            do {
                sb.append((char)ch);
            } while((ch = dis.read()) > 0);
            list.add(sb.toString());
        }

        headers = list.toArray(new String[list.size()]);
    }
    
    /**
     * Writes this header into the buffer.
     * 
     * @param buf the buffer to write the header into.
     */
    public void write(ByteBuffer buf) {
        buf.put(MAGIC.getBytes());
        buf.put(version);
        buf.putLong(size);
        buf.putLong(len);

        for (String header : headers) {
            buf.put(header.getBytes());
            buf.put((byte)0);
        }
        buf.put((byte)0);
    }
    
    public int findHeader(String header) {
        for (int i = 0, n = headers.length; i < n; i++) {
            if (header.equals(headers[i])) {
                return i;
            }
        }
        return -1;
    }
    
    public int getBlockHeaderLength() {
        return getBlockHeaderLength(headers);
    }

    public long getHeaderHash() {
        return getBlockHeaderHash(headers);
    }
    
    public static long getBlockHeaderHash(String[] headers) {
        long hash = 1125899906842597L;
        for (String header : headers) {
            for (int i = 0, n = header.length(); i < n; i++) {
                hash = (hash << 5) - hash + header.charAt(i);
            }
        }
        return hash;
    }
    
    public static int getBlockHeaderLength(String[] headers) {
        int len = 26; // 8 ('GecozBWT') + 1 (version) + 8 (size) + 8 (length) + 1 (last '\0') 
        for (String hdr : headers) {
            len += hdr.length() + 1; // '\0'
        }
        return len;
    }
}
