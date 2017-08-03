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

package es.elixir.bsc.ngs.nova.gzip;

import es.elixir.bsc.ngs.nova.io.DataReaderHelper;
import java.io.IOException;
import java.io.InputStream;
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
    
    public GZipHeader(InputStream in) throws IOException, ZipException {
        if (in.read() != ID1 ||
            in.read() != ID2) {
            throw new ZipException("invalid gzip header");
        }

        if (in.read() != CM) {
            throw new ZipException("unknown compression");
        }

        flg = in.read();
//        if (file.readUnsignedByte() != FLG) {
//            throw new ZipException("invalid extra flag for bgzf format");
//        }

        mtime = DataReaderHelper.readUnsignedInt(in); // !!!!
        xfl = in.read();
        os = in.read();

        if ((flg & FEXTRA) == 0) {
            xlen = 0;
            dsize = 0;
        } else {
            xlen = DataReaderHelper.readUnsignedShort(in);

            int bsize = 0;
            for (int len = xlen; len > 0;) {
                if (len < 4) {
                    throw new ZipException("invalid extra field descriptor");
                }

                final int si1 = in.read();
                final int si2 = in.read();
                final int slen = DataReaderHelper.readUnsignedShort(in);

                len -= 4;

                if (si1 == SI1 && si2 == SI2 && slen == SLEN) {
                    bsize = DataReaderHelper.readUnsignedShort(in);
                    break;
                }
            }

            dsize = bsize - xlen - 19;
        }
        
        if ((flg & FNAME) != 0) {
            StringBuilder str = new StringBuilder();
            int ch;
            while((ch = in.read()) > 0) {
                str.append((char)ch);
            }
        }
        
        if ((flg & FCOMMENT) != 0) {
            StringBuilder str = new StringBuilder();
            int ch;
            while((ch = in.read()) > 0) {
                str.append((char)ch);
            }
        }
        
        if ((flg & FHCRC) != 0) {
            DataReaderHelper.readUnsignedShort(in);
        }
    }
}
