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

package es.elixir.bsc.ngs.nova.sam.header;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

/**
 * @author Dmitry Repchevsky
 */

public abstract class AbstractHeaderLine {

    public abstract void write(PrintStream out) throws IOException;
    
    @Override
    public String toString() {
        try (ByteArrayOutputStream out = new ByteArrayOutputStream()) {
            write(new PrintStream(out, true, StandardCharsets.US_ASCII.name()));
            return out.toString(StandardCharsets.US_ASCII.name());
        } catch (IOException ex) {}
        return "";
    }
    
    /**
     * <p>
     * Parses the SAM header line looking for defined tags.
     * </p>
     * N.B. Method returns the same input tags array with found values or null 
     * if no value found for the tag.
     * 
     * @param line SAM header line to be parsed
     * @param tags SAM header tags to be found in the line
     * 
     * @return an array of requested tags with found values
     */
    public static String[] parseHeaderLine(final String line, final String[] tags) {
        
        final String[] elements = line.split("\t");
        if(elements.length == 0) {
            throw new IllegalArgumentException();
        }
        
        label:
        for (int i = 0, n = tags.length; i < n; i++) {
            for (int j = 0, m = elements.length; j < m; j++) {
                if (elements[j].length() > 3 && 
                    elements[j].startsWith(tags[i]) &&
                    elements[j].charAt(2) == ':') {
                    tags[i] = elements[j].substring(3);
                    continue label;
                }
            }
            tags[i] = null;
        }
        return tags;
    }

}
