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

package es.elixir.bsc.ngs.nova.io;

import java.io.IOException;

/**
 * 
 * @author Dmitry Repchevsky
 */

public interface BitInputStream {

    /**
     * Reads bits from the stream without moving the pointer.
     * 
     * @param nbits the number of bits to read
     * 
     * @return the long value with nbits bits. 
     * <b>NB:</b> the value is dirty! It may contain more bits than requested.
     * 
     * @throws IOException 
     */
    long peekBits(int nbits) throws IOException;
    
    /**
     * Discards bits from the stream.
     * 
     * @param nbits the number of bits to discard
     * 
     * @throws IOException 
     */
    void skipBits(int nbits) throws IOException;
    
    /**
     * Reads bits from the stream.
     * 
     * @param nbits the number of bits to read
     * 
     * @return the long value with nbits bits. 
     * <b>NB:</b> the value is dirty! It may contain more bits than requested.
     * 
     * @throws IOException 
     */
    long readBits(int nbits) throws IOException;
    
    /**
     * Aligns read position to the next byte.
     * 
     * @throws IOException 
     */
    void align() throws IOException;
}
