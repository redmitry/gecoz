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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
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
    
    public BAMFileReader(File fbam) throws IOException, DataFormatException {
        this(fbam, null);
    }

    public BAMFileReader(File fbam, File fbai) throws IOException, DataFormatException {
        
        if (fbai != null) {
            if (fbai.exists()) {
                try (FileInputStream in = new FileInputStream(fbai)) {
                    bai = new BAI(new BufferedInputStream(in));
                }
            } else {
                saveIndex(fbam, fbai);
            }
        } else {
            bai = makeIndex(fbam);
        }
        
        stream = new BAMFileInputStream(fbam);
    }
    
    public List<BAMAlignment> search(int idRef, int start, int end) throws IOException, DataFormatException {
        List<BAMAlignment> alignments = new ArrayList();
        
        Bin[] bins = bai.indexes[idRef];
        if (bins != null) {
            int[] nbins = BAI.reg2bins(start, end);
            for (int nbin : nbins) {
                final Bin bin = bins[nbin];
                final int n_chunk = bin.chunk_beg.length;
                for (int i = 0; i < n_chunk; i++) {
                    final long chunk_beg = bin.chunk_beg[i];
                    final long chunk_end = bin.chunk_end[i];

                    stream.move(chunk_beg);
                    while (stream.available() >= 0 && stream.index() < chunk_end) {
                        BAMAlignment a = BAMAlignment.decode(stream);
                        if (a.getPositionStart() < end && a.getPositionEnd() > start) {
                            alignments.add(a);
                        }
                    }
                }
            }
        }
        
        return alignments;
    }
    
    private BAI saveIndex(File fbam, File fbai) throws IOException, DataFormatException {
        bai = makeIndex(fbam);
        try (FileOutputStream out = new FileOutputStream(fbai)) {
            bai.save(out);
        }
        return bai;
    }
    
    private BAI makeIndex(File fbam) throws IOException, DataFormatException {
        try (BAMFileInputStream bam = new BAMFileInputStream(fbam)) {
            return new BAI(bam);
        }
    }
}
