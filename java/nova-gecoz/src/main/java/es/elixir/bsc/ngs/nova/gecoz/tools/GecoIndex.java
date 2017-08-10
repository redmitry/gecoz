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

package es.elixir.bsc.ngs.nova.gecoz.tools;

import es.elixir.bsc.ngs.nova.fasta.FastaFileReader;
import es.elixir.bsc.ngs.nova.fasta.FastaIterator;
import es.elixir.bsc.ngs.nova.fasta.FastaSequence;
import es.elixir.bsc.ngs.nova.gecoz.GecozFileWriter;
import es.elixir.bsc.ngs.nova.gecoz.GecozRefBlock;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Dmitry Repchevsky
 */
public class GecoIndex {
    
    static void index(Path ipath, Path opath, Path xpath, int sampling) {
        GecoIndex.index(ipath, opath, xpath, sampling, 1);
    }
    
    static void index(Path ipath, Path opath, Path xpath, int sampling, int threads) {

        Logger.getLogger(GecoIndex.class.getName()).log(Level.INFO, "analysing fasta file: {0} ...\n", ipath);

        final long t1 = System.nanoTime();

        TreeSet<GecozRefBlock> blocks = new TreeSet<>();
        
        FastaFileReader reader = new FastaFileReader(ipath, true);
        try (FastaIterator iter = reader.iterator()) {
            while(iter.hasNext()) {
                blocks.add(new GecozRefBlock(iter.next()));
            }
        } catch(Exception ex) {
            Logger.getLogger(GecoIndex.class.getName()).log(Level.SEVERE, "error reading file: {0}\n", ipath);
            Logger.getLogger(GecoIndex.class.getName()).log(Level.SEVERE, ex.getMessage());
            System.exit(1);
        }
        
        if (blocks.isEmpty()) {
            Logger.getLogger(GecoIndex.class.getName()).log(Level.SEVERE, "no data found in file: {0}\n", ipath);
            System.exit(1);
        }
        
        final int max_size = blocks.last().size();
        while (blocks.size() > 1) {
            GecozRefBlock first = blocks.pollFirst();
            GecozRefBlock second = blocks.pollFirst();
            final int size = first.size() + second.size();
            if (size > 0 && size <= max_size) {
                first.add(second.sequences);
                blocks.add(first);
            } else {
                blocks.add(first);
                blocks.add(second);
                break;
            }
        }
        
        // sort blocks to put those that have largest sequence first
        TreeSet<GecozRefBlock> sorted = new TreeSet(new Comparator<GecozRefBlock>() {
            @Override
            public int compare(GecozRefBlock o1, GecozRefBlock o2) {
                if (o1.sequences.first().length != o2.sequences.first().length) {
                    return o1.sequences.first().length > o2.sequences.first().length ? -1 : 1;
                }
                return o1.compareTo(o2);
            }
        });

        sorted.addAll(blocks);

        try (GecozFileWriter writer = new GecozFileWriter(opath, xpath, sampling, threads)) {
            for (GecozRefBlock block : sorted) {
                writeBlock(reader, writer, block);
            }
        } catch (Throwable th) {
            th.printStackTrace(System.err);
            System.exit(1);
        }

        final long t2 = System.nanoTime();
        Logger.getLogger(GecoIndex.class.getName()).log(Level.INFO, "finished in {0} ms.\n", ((t2 - t1)/1000000));
    }
    
    private static void writeBlock(FastaFileReader reader, GecozFileWriter writer, GecozRefBlock block) throws IOException {
        
        ByteBuffer buf;
        boolean attempt = true;
        do {
            try {
                buf = ByteBuffer.allocate(block.size());
                break;
            } catch(OutOfMemoryError th) {
                if (attempt = !attempt) {
                    throw th;                    
                }
                System.gc();
                Logger.getLogger(GecozFileWriter.class.getName()).log(Level.WARNING, "warning: low memory (free: {0} bytes)\n", Runtime.getRuntime().freeMemory());
            }
        } while(true);

        int i = 0;
        String[] headers = new String[block.sequences.size()];
        for (FastaSequence seq : block.sequences) {
            headers[i++] = seq.header;
            reader.read(buf, seq);
            buf.put((byte)0);
        }
        buf.rewind();

        writer.write(headers, buf);
    }
}
