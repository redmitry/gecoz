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

/**
 * @author Dmitry Repchevsky
 */

public class ReferenceLineAlternativeLocus {
    private final String chromosome;
    private final long start;
    private final long end;

    public ReferenceLineAlternativeLocus(
            final String chromosome, 
            final long start, 
            final long end) {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
    }

    public ReferenceLineAlternativeLocus(String substring) {
        String[] firstSplit = substring.split(":");
        if(firstSplit.length != 2){
            throw new IllegalArgumentException();
        }
        String[] secondSplit = firstSplit[1].split("-");
        if(secondSplit.length != 2) {
            throw new IllegalArgumentException();
        }
        chromosome = firstSplit[0];
        try {
            start = Long.parseLong(secondSplit[0]);
            end = Long.parseLong(secondSplit[1]);
        }catch (NumberFormatException e){
            throw new IllegalArgumentException(e);
        }
    }

    public String getChromosome() {
        return chromosome;
    }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return end;
    }

    @Override
    public String toString() {
        return chromosome + ':' + start + "-" + end;
    }
}
