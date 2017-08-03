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

package es.elixir.bsc.ngs.nova.bam;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import es.elixir.bsc.ngs.nova.io.DataReaderHelper;

/**
 * @author Dmitry Repchevsky
 */

public class Bin implements Comparable<Bin> {
        
    public final int bin;

    public long[] chunk_beg;
    public long[] chunk_end;

    protected Bin(int bin, long[] chunk_beg, long[] chunk_end) {
        this.bin = bin;
        this.chunk_beg = chunk_beg;
        this.chunk_end = chunk_end;
    }

    public Bin(InputStream in) throws IOException {
        bin = (int)DataReaderHelper.readUnsignedInt(in);

        final int n_chunk = (int)DataReaderHelper.readUnsignedInt(in);

        chunk_beg = new long[n_chunk];
        chunk_end = new long[n_chunk];

        for (int i = 0; i < n_chunk; i++) {
            chunk_beg[i] = DataReaderHelper.readLong(in);
            chunk_end[i] = DataReaderHelper.readLong(in);
        }
    }

    /**
     * Adds the chunk to the bin.
     * 
     * @param idx1 starting position of the chunk in the file
     * @param idx2 last position of the chunk in the file
     */
    public void add(final long idx1, final long idx2) {
        final int idx = chunk_beg.length;

        chunk_beg = Arrays.copyOf(chunk_beg, idx + 1);
        chunk_end = Arrays.copyOf(chunk_end, idx + 1);

        chunk_beg[idx] = idx1;
        chunk_end[idx] = idx2;
    }

    /**
     * Merges the chunk with previous one if within the same gzip chunk.
     * Adds otherwise.
     * 
     * @param idx1 starting position of the chunk in the file
     * @param idx2 last position of the chunk in the file
     */
    public void merge(final long idx1, final long idx2) {
        final int idx = chunk_end.length - 1;

        if (chunk_end[idx] >= idx1) {
            if (chunk_end[idx] < idx2) {
                chunk_end[idx] = idx2;
            }
        } else if (chunk_end[idx] >> 16 == idx2 >> 16) {
            chunk_end[idx] = idx2;
        } else {
            add(idx1, idx2);
        }            
    }

    @Override
    public int compareTo(Bin o) {
        return this.bin == o.bin ? 0 : this.bin > o.bin ? 1 : -1;
    }
}