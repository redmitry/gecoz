/**
 * *****************************************************************************
 * Copyright (C) 2016 Spanish National Bioinformatics Institute (INB) and
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

package es.elixir.bsc.ngs.nova.bam;

/**
 * @author Dmitry Repchevsky
 */

public class Quality {
    private final byte[] qual;
    
    public Quality(String quality) {
        qual = Quality.parse(quality);
    }

    public Quality(byte[] quality) {
        qual = quality;
    }
    
    @Override
    public String toString() {
        return qual != null ? parse(qual) : "";
    }
    
    public static String parse(byte[] quality) {
        StringBuilder sb = new StringBuilder(quality.length);
        for (int i = 0, n = quality.length; i < n; i++) {
            sb.append((char)(quality[i] + 33));
        }
        return sb.toString();
    }
    
    public static byte[] parse(String quality) {
        final byte[] qual = new byte[quality.length()];
        for (int i = 0, n = quality.length(); i < n; i++) {
            qual[i] = (byte) (quality.charAt(i) - 33);
        }
        return qual;
    }
}
