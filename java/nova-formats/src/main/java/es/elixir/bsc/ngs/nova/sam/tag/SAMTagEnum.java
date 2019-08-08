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

package es.elixir.bsc.ngs.nova.sam.tag;

import es.elixir.bsc.ngs.nova.sam.SAMTag;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * @author Dmitry Repchevsky
 */

public enum SAMTagEnum {
    
    AM('i'), AS('i'), BC('Z'), BQ('Z'), BZ('Z'), CB('Z'), CC('Z'), CG('B'),
    CM('i'), CO('Z'), CP('i'), CQ('Z'), CR('Z'), CS('Z'), CT('Z'), CY('Z'),
    E2('Z'), FI('i'), FS('Z'), FZ('B'), GC('?'), GQ('?'), GS('?'), H0('i'), 
    H1('i'), H2('i'), HI('i'), IH('i'), LB('Z'), MC('Z'), MD('Z'), MF('?'),
    MI('Z'), MQ('i'), NH('i'), NM('i'), OA('Z'), OC('Z'), OP('i'), OQ('Z'),
    OX('Z'), PG('Z'), PQ('i'), PT('Z'), PU('Z'), Q2('Z'), QT('Z'), QX('Z'),
    R2('Z'), RG('Z'), RT('?'), RX('Z'), S2('?'), SA('Z'), SM('i'), SQ('?'),
    TC('i'), U2('Z'), UQ('i'), ;

    public final char type;
    
    SAMTagEnum(final char type) {
        this.type = type;
    }
    
    public static <T extends SAMTag> T decode(final SAMTagEnum tag, final ByteBuffer buf) {
        return (T)decode(tag.name(), buf);
    }

    public static SAMTag decode(final String tag, final ByteBuffer buf) {
        final char val_type = (char)buf.get();
        switch(val_type) {
            case 'c' : 
            case 'C' : return decode(tag, val_type, buf.get());
            case 's' : 
            case 'S' : return decode(tag, val_type, buf.getShort());
            case 'i' : 
            case 'I' : return decode(tag, val_type, buf.getInt());
            case 'Z' : byte ch0;
                       final StringBuilder sb1 = new StringBuilder();
                       while ((ch0 = buf.get()) != 0) {
                           sb1.append((char)ch0);
                       }
                       return decode(tag, val_type, sb1.toString());
            case 'B' : final byte a_type = buf.get();
                       final int count = buf.getInt();
                       switch (a_type) {
                           case 'c':
                           case 'C': final byte[] b = new byte[count];
                                     buf.get(b);
                                     return new AnySAMTag(tag, val_type, b);
                           case 's': 
                           case 'S': final short[] s = new short[count];
                                     for (int c = 0; c < count; c++) {
                                         s[c] = buf.getShort();
                                     }
                                     return new AnySAMTag(tag, val_type, s);
                           case 'i':
                           case 'I': final int[] ui = new int[count];
                                     for (int c = 0; c < count; c++) {
                                         ui[c] = buf.getInt();
                                     }
                                     return decode(tag, val_type, ui);
                       }
        }
        
        return null;
    }

    private static SAMTag decode(final String tag, final char val_type, final int b) {
        switch(tag) {
            case "AM" : return new AM(b);
            case "AS" : return new AS(b);
            case "CM" : return new CM(b);
            case "CP" : return new CP(b);
            case "FI" : return new FI(b);
            case "H0" : return new H0(b);
            case "H1" : return new H1(b);
            case "H2" : return new H2(b);
            case "HI" : return new HI(b);
            case "IH" : return new IH(b);
            case "MQ" : return new MQ(b);
            case "NH" : return new NH(b);
            case "NM" : return new NM(b);
            case "OP" : return new OP(b);
            case "PQ" : return new PQ(b);
            case "SM" : return new SM(b);
        }
        return new AnySAMTag(tag, val_type, b);
    }
    
    private static SAMTag decode(final String tag, final char val_type, final String z) {
        switch(tag) {
            case "BC" : return new BC(z);
            case "BQ" : return new BQ(z);
            case "BZ" : return new BZ(z);
            case "CB" : return new CB(z);
            case "CC" : return new CC(z);
            case "CO" : return new CO(z);
            case "CQ" : return new CQ(z);
            case "CR" : return new CR(z);
            case "CS" : return new CS(z);
            case "CT" : return new CT(z);
            case "CY" : return new CY(z);
            case "E2" : return new E2(z);
            case "FS" : return new FS(z);
            case "LB" : return new LB(z);
            case "MC" : return new MC(z);
            case "MD" : return new MD(z);
            case "MI" : return new MI(z);
            case "OA" : return new OA(z);
            case "OC" : return new OC(z);
            case "OQ" : return new OQ(z);
            case "OX" : return new OX(z);
            case "PG" : return new PG(z);
            case "PT" : return new PT(z);
            case "PU" : return new PU(z);
            case "Q2" : return new Q2(z);
            case "QT" : return new QT(z);
            case "QX" : return new QX(z);
            case "R2" : return new R2(z);
            case "RG" : return new RG(z);
            case "RX" : return new RX(z);
            case "SA" : return new SA(z);
            case "U2" : return new U2(z);
        }
        return new AnySAMTag(tag, val_type, z);
    }

    private static SAMTag decode(final String tag, final char val_type, final int[] i) {
        switch(tag) {
            case "CG" : return new CG(i);
            case "FZ" : return new FZ(i);
            
        }
        return new AnySAMTag(tag, val_type, i);
    }

    public static byte[] encode(final SAMTag tag) {
        final String name = tag.getTagName();
        
        byte[] data;
        
        final char type = tag.getTagType();
        final Object value = tag.getTagValue();
        if (value.getClass().isArray()) {
            switch(type) {
                case 'c':
                case 'C': final byte[] b = (byte[])value;
                          final ByteBuffer bbuf = ByteBuffer.wrap(data = new byte[8 + b.length]).order(ByteOrder.LITTLE_ENDIAN);
                          bbuf.position(8);
                          bbuf.put(b);
                          break; 
                case 's':
                case 'S': final short[] s = (short[])value;
                          final ByteBuffer sbuf = ByteBuffer.wrap(data = new byte[8 + s.length * 2]).order(ByteOrder.LITTLE_ENDIAN);
                          sbuf.position(8);
                          sbuf.asShortBuffer().put(s);
                          break;
                case 'i':
                case 'I': final int[] i = (int[])value;
                          final ByteBuffer ibuf = ByteBuffer.wrap(data = new byte[8 + i.length * 4]).order(ByteOrder.LITTLE_ENDIAN);
                          ibuf.position(8);
                          ibuf.asIntBuffer().put(i);
                          break;
                default: return null;
            }
        } else {

            switch(type) {
                case 'c':
                case 'C': ByteBuffer.wrap(data = new byte[4]).order(ByteOrder.LITTLE_ENDIAN).put(3, (byte)value);
                          break;
                case 's': 
                case 'S': ByteBuffer.wrap(data = new byte[5]).order(ByteOrder.LITTLE_ENDIAN).putShort(3, (short)value);
                          break;
                case 'i': 
                case 'I': ByteBuffer.wrap(data = new byte[7]).order(ByteOrder.LITTLE_ENDIAN).putInt(3, (int)value);
                          break;
                case 'Z': final ByteBuffer buf = ByteBuffer.wrap(data = new byte[4 + value.toString().length()]).order(ByteOrder.LITTLE_ENDIAN);
                          buf.position(3);
                          buf.put(value.toString().getBytes());
                          break;
                case 'f': ByteBuffer.wrap(data = new byte[7]).order(ByteOrder.LITTLE_ENDIAN).putFloat(3, (float)value);
                          break;
                          

                default: return null;
            }
            data[2] = (byte)type;
        }
        
        data[0] = (byte)name.charAt(0);
        data[1] = (byte)name.charAt(1);
        
        return data;
    }
    
    public byte[] encode(final long b) {
        switch(type) {
            case 'c' : 
            case 'C' : return new byte[] 
                        {(byte)name().charAt(0), (byte)name().charAt(1), (byte)type, 
                         (byte)b};
            case 's' : 
            case 'S' : return new byte[] 
                        {(byte)name().charAt(0), (byte)name().charAt(1), (byte)type, 
                         (byte)(b & 0xFF), (byte)(b >> 8)};
            case 'i' : 
            case 'I' : return new byte[] 
                        {(byte)name().charAt(0), (byte)name().charAt(1), (byte)type, 
                         (byte)(b & 0xFF), (byte)(b >> 8), (byte)(b >> 16), (byte)(b >> 32)};
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
