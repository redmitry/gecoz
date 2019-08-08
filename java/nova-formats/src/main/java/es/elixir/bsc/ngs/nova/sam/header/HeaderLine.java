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

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * @author Dmitry Repchevsky
 */

public class HeaderLine extends AbstractHeaderLine {
    
    public final static String TAG = "@HD";

    public final static String[] TAGS = {"VN", "SO", "GO", "SS"};
    
    public final String version;
    public final SortingOrder sortingOrder;
    public final AlignmentsGrouping alignmentsGrouping;
    public final String subSorting;

    public HeaderLine(
            final String version,
            final SortingOrder sortingOrder,
            final AlignmentsGrouping alignmentsGrouping,
            final String subSorting) {
        this.version = version;
        this.sortingOrder = sortingOrder;
        this.alignmentsGrouping = alignmentsGrouping;
        this.subSorting = subSorting;
    }

    public HeaderLine(final String line) {

        final String[] tags = parseHeaderLine(line, Arrays.copyOf(TAGS, TAGS.length));

        version = tags[0];
        sortingOrder = tags[1] == null ? null : SortingOrder.valueOf(tags[1].toUpperCase());
        alignmentsGrouping = tags[2] == null ? null : AlignmentsGrouping.valueOf(tags[2].toUpperCase());
        subSorting = tags[3];
    }
    
    @Override
    public void write(final PrintStream out) throws IOException {
        out.append(TAG);
        
        if (version != null) {
            out.append("\tVN:").append(version);
        }
        if(sortingOrder != null) {
            out.append("\tSO:").print(sortingOrder);
        }
        if(alignmentsGrouping != null) {
            out.append("\tGO:").print(alignmentsGrouping);
        }
        if(subSorting != null) {
            out.append("\tSS:").append(subSorting);
        }
    }
}
