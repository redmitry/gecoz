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

package es.elixir.bsc.ngs.nova.gzip;

import es.elixir.bsc.ngs.nova.io.DataReaderHelper;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.util.zip.ZipException;

/**
 * @author Dmitry Repchevsky
 */

public class GZipHeader {
    public final static int ID1 = 0x1F;
    public final static int ID2 = 0x8B;
    public final static int CM = 0x08;
    
    public final static int FHCRC = 0b10;
    public final static int FEXTRA = 0b100;
    public final static int FNAME = 0b1000;
    public final static int FCOMMENT = 0b10000;
    
    public final int flg;
    public final long mtime;
    public final int xfl;
    public final int os;
    public final int xlen;

    public final static int SI1 = 0x42;
    public final static int SI2 = 0x43;
    public final static int SLEN = 0x02;
    public final int dsize;
    
    public final String file_name;
    public final String file_comment;
    
    public GZipHeader() {
        flg = 0;
        mtime = 0; // no timestamp;
        xfl = 2;   // maximum compression
        os = 255;  // unknown
        xlen = 0;
        dsize = 0;
        file_name = null;
        file_comment = null;
    }

    public GZipHeader(final long mtime,
                      final String file_name,
                      final String file_comment,
                      final int bsize) {
        this.mtime = mtime; // no timestamp;
        xfl = 2;   // maximum compression
        os = 255;  // unknown
        this.file_name = file_name;
        this.file_comment = file_comment;
        
        flg = (file_name == null ? 0 : FNAME) | (file_comment == null ? 0 : FCOMMENT) | FEXTRA;
        xlen = 6;
        int sz = 0;
        if (bsize > 25) {
            if (file_name != null) {
                sz = sz - file_name.length() - 1;
            }
            if (file_comment != null) {
                sz = sz - file_comment.length() - 1;
            }
            sz -= 25; // BSIZE-XLEN-19
        }
        dsize = sz > 0 ? sz : 0;
    }

    public GZipHeader(final InputStream in) throws IOException {
        final int id1 = in.read();
        if (id1 < 0) {
            throw new EOFException();
        }
        if (id1 != ID1 ||
            in.read() != ID2) {
            throw new ZipException("invalid gzip header");
        }

        if (in.read() != CM) {
            throw new ZipException("unknown compression");
        }

        flg = in.read();

        mtime = DataReaderHelper.readUnsignedInt(in); // !!!!
        xfl = in.read();
        os = in.read();

        int bsize = 0;
        if ((flg & FEXTRA) == 0) {
            xlen = 0;
        } else {
            xlen = DataReaderHelper.readUnsignedShort(in);
            for (int len = xlen; len > 0;) {
                if (len < 4) {
                    throw new ZipException("invalid extra field descriptor");
                }

                final int si1 = in.read();
                final int si2 = in.read();
                final int slen = DataReaderHelper.readUnsignedShort(in);

                len -= 4;

                if (si1 == SI1 && si2 == SI2 && slen == SLEN) {
                    bsize = DataReaderHelper.readUnsignedShort(in); // Block SIZE minus 1
                    break;
                }
            }
            bsize -= xlen;
        }
        
        if ((flg & FNAME) != 0) {
            final StringBuilder str = new StringBuilder();
            int ch;
            while((ch = in.read()) > 0) {
                str.append((char)ch);
            }
            file_name = str.toString();
            bsize = bsize - file_name.length() - 1;
        } else {
            file_name = null;
        }
        
        if ((flg & FCOMMENT) != 0) {
            final StringBuilder str = new StringBuilder();
            int ch;
            while((ch = in.read()) > 0) {
                str.append((char)ch);
            }
            file_comment = str.toString();
            bsize = bsize - file_comment.length() - 1;
        } else {
            file_comment = null;
        }
        
        if ((flg & FHCRC) != 0) {
            DataReaderHelper.readUnsignedShort(in);
            bsize -= 2;
        }
        
        // 18 == ID1 + ID2 + CM + FLG + MTIME + XFL + OS + XLEN + CRC32 + ISIZE
        if (bsize <= 19) {
            dsize = 0; // BSIZE is not set or wrong
        } else {
            dsize = bsize - 19; // BSIZE-XLEN-19
        }
    }
    
    public void write(final OutputStream out) throws IOException {
        out.write(ID1);
        out.write(ID2);
        out.write(CM);
        out.write(flg);
        
        // write time
        out.write((int)(mtime & 0xFF));
        out.write((int)((mtime >>> 8) & 0xFF));
        out.write((int)((mtime >>> 16) & 0xFF));
        out.write((int)((mtime >>> 24) & 0xFF));
        
        out.write(xfl);
        out.write(os);
        
        int bsize = dsize + xlen + 19; // DSIZE+XLEN+19
        if (file_name != null) {
            bsize = bsize + file_name.length() + 1;
        }
        if (file_comment != null) {
            bsize = bsize + file_comment.length() + 1;
        }

        if ((flg & FEXTRA) != 0) {
            // write xlen
            out.write(xlen & 0xFF);
            out.write((xlen >>> 8) & 0xFF);
            
            out.write(SI1);
            out.write(SI2);
            
            out.write(SLEN);
            out.write(0);

            if (dsize == 0) {
                out.write(0);
                out.write(0);
            } else {
                out.write(bsize & 0xFF);
                out.write((bsize >>> 8) & 0xFF);
            }
        }
        
        if (file_name != null) {
            out.write(file_name.getBytes(StandardCharsets.ISO_8859_1));
            out.write(0);
        }
        
        if (file_comment != null) {
            out.write(file_comment.getBytes(StandardCharsets.ISO_8859_1));
            out.write(0);
        }
    }
    
    public int getHeaderSize() {
        int bsize = xlen + 18;
        if (file_name != null) {
            bsize = bsize + file_name.length() + 1;
        }
        if (file_comment != null) {
            bsize = bsize + file_comment.length() + 1;
        }
        return bsize;
    }
}
