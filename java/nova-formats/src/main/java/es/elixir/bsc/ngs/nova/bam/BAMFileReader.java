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
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.DataFormatException;

/**
 * The main class to work with a BAM file.
 * 
 * @author Dmitry Repchevsky
 */

public class BAMFileReader {
    
    private BAI bai;
    private BAMFileInputStream stream;
    
    public BAMFileReader(final Path fbam) throws IOException, DataFormatException {
        this(fbam, null);
    }

    public BAMFileReader(final Path fbam, final Path fbai) throws IOException, DataFormatException {
        
        if (fbai != null) {
            if (Files.exists(fbai) && Files.isRegularFile(fbai)) {
                try (InputStream in = Files.newInputStream(fbai)) {
                    bai = new BAI(in);
                }
            } else {
                saveIndex(fbam, fbai);
            }
        } else {
            bai = makeIndex(fbam);
        }
        
        stream = new BAMFileInputStream(fbam);
    }
    
    public List<BAMAlignment> search(final int idRef, 
                                     final int start, 
                                     final int end) throws IOException, DataFormatException {

        final List<BAMAlignment> alignments = new ArrayList();
        
        final Bin[] bins = bai.indexes[idRef];
        if (bins != null) {
            final int[] nbins = BAI.reg2bins(start, end);
            for (int i = 0; i < nbins.length; i++) {
                final Bin bin = bins[nbins[i]];
                final int n_chunk = bin.chunk_beg.length;
                for (int j = 0; j < n_chunk; j++) {
                    final long chunk_beg = bin.chunk_beg[j];
                    final long chunk_end = bin.chunk_end[j];

                    stream.move(chunk_beg);
                    while (stream.available() >= 0 && stream.index() < chunk_end) {
                        final BAMAlignment a = BAMAlignment.decode(stream);
                        if (a.getPositionStart() < end && a.getPositionEnd() > start) {
                            alignments.add(a);
                        }
                    }
                }
            }
        }
        
        return alignments;
    }
    
    private BAI saveIndex(final Path fbam, final Path fbai) throws IOException, DataFormatException {
        bai = makeIndex(fbam);
        try (OutputStream out = Files.newOutputStream(fbai)) {
            bai.save(out);
        }
        return bai;
    }
    
    private BAI makeIndex(final Path fbam) throws IOException, DataFormatException {
        try (BAMFileInputStream bam = new BAMFileInputStream(fbam)) {
            return new BAI(bam);
        }
    }
}
