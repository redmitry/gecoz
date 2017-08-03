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

package es.elixir.bsc.ngs.nova.sam;

import java.nio.ByteBuffer;
import javax.xml.bind.DatatypeConverter;

/**
 * @author Dmitry Repchevsky
 */

public enum SAMTag {
    
    AM('i'), AS('i'), BC('Z'), BQ('Z'), CC('Z'), CM('i'), CO('Z'), CP('i'),
    CQ('Z'), CS('Z'), CT('Z'), E2('Z'), FI('i'), FS('Z'), FZ('B'), LB('Z'),
    H0('i'), H1('i'), H2('i'), HI('i'), IH('i'), MC('Z'), MD('Z'), MQ('i'),
    NH('i'), NM('i');

    public final byte type;
    
    SAMTag(final char type) {
        this.type = (byte)type;
    }
    
    public static Object decode(ByteBuffer buf) {
        final char val_type = (char)buf.get();
        switch(val_type) {
            case 'A' : return Character.toString((char)buf.get());
            case 'c' : return buf.get();
            case 'C' : return buf.get() & 0xFF;
            case 's' : return buf.getShort();
            case 'S' : return buf.getShort() & 0xFFFF;
            case 'i' : return buf.getInt();
            case 'I' : return buf.getInt() & 0xFFFFFFFFL;
            case 'f' : return buf.getFloat();
            case 'Z' : byte ch0;
                       StringBuilder sb1 = new StringBuilder();
                       while ((ch0 = buf.get()) != 0) {
                           sb1.append((char)ch0);
                       }
                       return sb1.toString();
            case 'H' : byte ch1;
                       StringBuilder sb2 = new StringBuilder();
                       while ((ch1 = buf.get()) != 0) {
                           sb2.append((char)ch1);
                       }
                       return DatatypeConverter.parseHexBinary(sb2.toString());
            case 'B' : final byte a_type = buf.get();
                       final int len = buf.getInt();
                       String[] val = new String[len];
                       for (int i = 0; i < len; i++) {
                           val[i] = decode(buf, a_type);
                       }
        }
        
        return null;
    }
    
    private static String decode(final ByteBuffer buf, final byte val_type) {
        switch(val_type) {
            case 'c' : return Byte.toString(buf.get());
            case 'C' : return Integer.toString(buf.get() & 0xFF);
            case 's' : return Short.toString(buf.getShort());
            case 'S' : return Integer.toString(buf.getShort() & 0xFFFF);
            case 'i' : return Integer.toString(buf.getInt());
            case 'I' : return Long.toString(buf.getInt() & 0xFFFFFFFFL);
            case 'f' : return Float.toString(buf.getFloat());
        }
        return "";
    }
}
