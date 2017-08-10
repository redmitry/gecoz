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

package es.elixir.bsc.ngs.nova.fasta;

/**
 * Truncated Fasta Sequence object that has no data
 * @author Dmitry Repchevsky
 */

public class TFastaSequence implements Comparable<TFastaSequence> {

    public final String header;
    public final int length;
    public final boolean multiline;
    
    public TFastaSequence(final String header, final int length, final boolean multiline) {
        this.header = header;
        this.length = length;
        this.multiline = multiline;
    }
    
    @Override
    public int compareTo(TFastaSequence o) {
        if (this.length != o.length) {
            return this.length > o.length ? -1 : 1;
        }

        return this.header.compareTo(o.header);
    }
}
