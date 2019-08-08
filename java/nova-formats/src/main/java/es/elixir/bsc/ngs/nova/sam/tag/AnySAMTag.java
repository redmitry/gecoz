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

/**
 * @author Dmitry Repchevsky
 */

public class AnySAMTag implements SAMTag {

    public final String name;
    public final char type;
    
    /**
     * The tag value which type depends on the TAG {@link #type}:
     * 
     * <pre>
     * 'c'  - Byte
     * 'C'  - Byte (0-255)
     * 's'  - Short
     * 'S'  - Short (0 - 65535)
     * 'i'  - Integer
     * 'I'  - Integer (0 - 2^32-1)
     * 'f'  - Float
     * 'A'  - Character
     * 'Z'  - String
     * 'H'  - String
     * 'Bc' - byte[]
     * 'BC' - byte[]
     * 'Bs' - short[]
     * 'BS' - short[]
     * 'Bi' - int[]
     * 'BI' - int[]
     * </pre>
     */
    public final Object value;
        
    public AnySAMTag(final String name, final char type, final int number) {
        if (name == null || name.length() != 2) {
            throw new IllegalArgumentException();
        }
        this.name = name;
        
        switch(type) {
            case 'c':
            case 'C': this.value = (byte)number; break;
            case 's':
            case 'S': this.value = number & 0xFFFF; break;
            case 'i':
            case 'I': this.value = (int)number; break;
            default: throw new IllegalArgumentException();
        }
        
        this.type = type;
    }

    public AnySAMTag(final String name, final char type, float decimal) {
        if (name == null || name.length() != 2) {
            throw new IllegalArgumentException();
        }
        this.name = name;
        
        switch(type) {
            case 'f': this.value = decimal; break;
            default: throw new IllegalArgumentException();
        }
        
        this.type = type;
    }

    public AnySAMTag(final String name, final char type, char character) {
        if (name == null || name.length() != 2) {
            throw new IllegalArgumentException();
        }
        this.name = name;
        
        switch(type) {
            case 'A': this.value = character; break;
            default: throw new IllegalArgumentException();
        }
        
        this.type = type;
    }
    
    public AnySAMTag(final String name, final char type, final String text) {
        if (name == null || name.length() != 2) {
            throw new IllegalArgumentException();
        }
        this.name = name;
        
        switch(type) {
            case 'Z': this.value = text; break;
            default: throw new IllegalArgumentException();
        }
        
        this.type = type;
    }

    public AnySAMTag(final String name, final char type, final byte[] array) {
        if (name == null || name.length() != 2) {
            throw new IllegalArgumentException();
        }
        this.name = name;
        
        switch(type) {
            case 'c':
            case 'C': this.value = array;
                      break;
            case 'H': this.value = new String(array);
                      break;
            default: throw new IllegalArgumentException();
        }

        this.type = type;
    }

    public AnySAMTag(final String name, final char type, final short[] array) {
        if (name == null || name.length() != 2) {
            throw new IllegalArgumentException();
        }
        this.name = name;
        
        switch(type) {
            case 's':
            case 'S': this.value = array;
                      break;
            default: throw new IllegalArgumentException();
        }

        this.type = type;
    }

    public AnySAMTag(final String name, final char type, final int[] array) {
        if (name == null || name.length() != 2) {
            throw new IllegalArgumentException();
        }
        this.name = name;
        
        switch(type) {
            case 'i':
            case 'I': this.value = array;
                      break;
            default: throw new IllegalArgumentException();
        }
        
        this.type = type;
    }

    @Override
    public char getTagType() {
        return type;
    }

    @Override
    public String getTagName() {
        return name;
    }

    @Override
    public Object getTagValue() {
        return value;
    }
}
