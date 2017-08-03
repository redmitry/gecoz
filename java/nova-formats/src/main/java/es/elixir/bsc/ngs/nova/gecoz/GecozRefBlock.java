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

package es.elixir.bsc.ngs.nova.gecoz;

import es.elixir.bsc.ngs.nova.fasta.FastaSequence;
import java.util.Set;
import java.util.TreeSet;

/**
 * Auxiliary class to fuse short sequences in a block.
 * 
 * @author Dmitry Repchevsky
 */

public class GecozRefBlock implements Comparable<GecozRefBlock>{
    
    private int size;
    public final TreeSet<FastaSequence> sequences = new TreeSet<>();
    
    public GecozRefBlock(FastaSequence sequence) {
        sequences.add(sequence);
        size = sequence.length + 1;
    }

    public void add(FastaSequence sequence) {
        sequences.add(sequence);
        size += sequence.length + 1;
    }

    public void add(Set<FastaSequence> sequences) {
        for (FastaSequence seq : sequences) {
            add(seq);
        }
    }
    
    public int size() {
        return size;
    }

    @Override
    public int compareTo(GecozRefBlock o) {
        if (this.size != o.size) {
            return this.size > o.size ? 1 : -1;
        }

        return this.sequences.first().compareTo(o.sequences.first());
    }
}
