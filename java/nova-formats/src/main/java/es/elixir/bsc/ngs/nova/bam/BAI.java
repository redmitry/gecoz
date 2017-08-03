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
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import java.util.zip.DataFormatException;
import es.elixir.bsc.ngs.nova.io.DataReaderHelper;

/**
 * @author Dmitry Repchevsky
 */

public class BAI {

    public final Bin[][] indexes;
    public final long[][] offsets;

    /**
     * Reads the BAI index from the input stream.
     * 
     * @param in the input stream with BAI index
     * 
     * @throws IOException
     * @throws DataFormatException 
     */
    public BAI(InputStream in) throws IOException, DataFormatException {
        if (in.read() != 'B' ||
            in.read() != 'A' ||
            in.read() != 'I' ||
            in.read() != 01) {
            throw new DataFormatException("invalid BAI header");
        }

        final int n_ref = (int) DataReaderHelper.readUnsignedInt(in);
        
        final TreeMap<Integer,Bin> map[] = new TreeMap[n_ref]; // bins;

        offsets = new long[n_ref][];
       
        for (int i = 0; i < n_ref; i++) {
            final int n_bin = (int) DataReaderHelper.readUnsignedInt(in);
            if (n_bin > 0) {
                map[i] = new TreeMap();
                for (int j = 0; j < n_bin; j++) {
                    Bin bin = new Bin(in);
                    map[i].put(bin.bin, bin);
                }
            }
            final int n_intv = (int) DataReaderHelper.readUnsignedInt(in);
            if (n_intv > 0) {
                offsets[i] = new long[n_intv];
                for (int j = 0; j < n_intv; j++) {
                    offsets[i][j] = DataReaderHelper.readLong(in);
                }
            }
        }
        
        indexes = indexes(map);
    }
    
    /**
     * Creates a new BAM index from the BAM file stream.
     * The stream must be already positioned to the BAM records 
     * to start reading (the BAM header must be already read).
     * 
     * @param in the BAM file to create the index
     * 
     * @throws IOException
     * @throws DataFormatException 
     */
    public BAI(BAMFileInputStream in) throws IOException, DataFormatException {
        final int nRef = in.getRefCount();
        
        final TreeMap<Integer,Bin> map[] = new TreeMap[nRef]; // bins;
        
        final BAMAlignment linear[] = new BAMAlignment[nRef];
        final List<Long> list[] = new ArrayList[nRef];

        while(in.available() >= 0) {
            
            final long chunk_beg = in.index();
            BAMAlignment a = BAMAlignment.decode(in);
            final long chunk_end = in.index();
            
            if (a.refID >= nRef) {
                continue;
            }
            
            final int ref = a.refID;
            
            final int beg = a.getPositionStart();
            final int end = a.getPositionEnd();

            final int nbin = reg2bin(beg - 1, end > 0 ? end : beg);

            if (map[ref] == null) {
                Bin bin = new Bin(nbin, new long[] {chunk_beg}, new long[] {chunk_end});
                map[ref] = new TreeMap();
                map[ref].put(nbin, bin);
            } else {
                Bin bin = map[ref].get(nbin);
                if (bin == null) {
                    bin = new Bin(nbin, new long[] {chunk_beg}, new long[] {chunk_end});
                    map[ref].put(nbin, bin);
                } else {
                    bin.merge(chunk_beg, chunk_end);
                }
            }

            final int lsegment = beg >> 14; // 16384
            final int rsegment = end >> 14;

            if (list[ref] == null) {
                list[ref] = new ArrayList();
                for (int i = 0; i < lsegment; i++) {
                    list[ref].add(0L);
                }
                for (int i = lsegment; i <= rsegment; i++) {
                    list[ref].add(chunk_beg);
                }
                linear[ref] = a;
            } else if (list[ref].size() <= rsegment) {
                for (int i = list[ref].size(); i <= lsegment; i++) {
                    list[ref].add(chunk_beg);
                }
                for (int i = lsegment; i < rsegment; i++) {
                    list[ref].add(chunk_beg);
                }
                linear[ref] = a;
            } else if (end > linear[ref].getPositionEnd()) {
                linear[ref] = a;
            }
        }
        
        offsets = new long[list.length][];
        for (int i = 0, n = offsets.length; i < n; i++) {
            if (list[i] != null) {
                offsets[i] = new long[list[i].size()];
                for (int j = 0, m = offsets[i].length; j < m; j++) {
                    offsets[i][j] = list[i].get(j);
                }
            }
        }

        indexes = indexes(map);
    }

    public void save(OutputStream out) throws IOException {
        out.write('B');
        out.write('A');
        out.write('I');
        out.write(01);
        
        ByteBuffer u32 = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        ByteBuffer u64 = ByteBuffer.allocate(8).order(ByteOrder.LITTLE_ENDIAN);

        final int n_ref = offsets.length;
        
        out.write(u32.putInt(0, n_ref).array());
        
        for (int i = 0; i < n_ref; i++) {
            Bin[] bins = indexes[i];
            if (bins == null) {
                out.write(u32.putInt(0, 0).array());
                out.write(u32.putInt(0, 0).array());
            } else {
                List<Bin> list = new LinkedList(Arrays.asList(bins));
                list.removeAll(Collections.singleton(null));
                out.write(u32.putInt(0, list.size()).array());
                for (Bin bin : list) {
                    out.write(u32.putInt(0, bin.bin).array());
                    final int n_chunk = bin.chunk_beg.length;
                    out.write(u32.putInt(0, n_chunk).array()); 
                    for (int j = 0; j < n_chunk; j++) {
                        out.write(u64.putLong(0, bin.chunk_beg[j]).array()); // 
                        out.write(u64.putLong(0, bin.chunk_end[j]).array());
                    }
                }
                final int i_intv = offsets[i].length;
                out.write(u32.putInt(0, i_intv).array());
                for (int j = 0; j < i_intv; j++) {
                    out.write(u64.putLong(0, offsets[i][j]).array()); 
                }
            }
        }
    }

    private Bin[][] indexes(TreeMap<Integer,Bin> map[]) {
        Bin[][] arr = new Bin[map.length][];
        for (int i = 0, n = arr.length; i < n; i++) {
            if (map[i] != null) {
                arr[i] = new Bin[map[i].lastEntry().getKey() + 1];
                Iterator<Bin> iter = map[i].values().iterator();
                while(iter.hasNext()) {
                    Bin bin = iter.next();
                    arr[i][bin.bin] = bin;
                }
            }
        }
        return arr;
    }

    public static int reg2bin(int start, int end) {
        end--;
        if (start >> 14 == end >> 14) {
            return ((1 << 15) - 1) / 7 + (start >> 14);
        }
        if (start >> 17 == end >> 17) {
            return ((1 << 12) - 1) / 7 + (start >> 17);
        }
        if (start >> 20 == end >> 20) {
            return ((1 << 9) - 1) / 7 + (start >> 20);
        }
        if (start >> 23 == end >> 23) {
            return ((1 << 6) - 1) / 7 + (start >> 23);
        }
        if (start >> 26 == end >> 26) {
            return ((1 << 3) - 1) / 7 + (start >> 26);
        }
        return 0;
    }
    
    public static int[] reg2bins(int start, int end) {
        int idx = 0;
        int[] bins = new int[1];
        for (int i = (start >> 26) + 1, n = (end >> 26) + 1; i <= n; i++) {
            bins = addBin(bins, idx++, i);
        }
        for (int i = (start >> 23) + 9, n = (end >> 23) + 9; i <= n; i++) {
            bins = addBin(bins, idx++, i);
        }
        for (int i = (start >> 20) + 73, n = (end >> 20) + 73; i <= n; i++) {
            bins = addBin(bins, idx++, i);
        }
        for (int i = (start >> 17) + 585, n = (end >> 17) + 585; i <= n; i++) {
            bins = addBin(bins, idx++, i);
        }
        for (int i = (start >> 14) + 4681, n = (end >> 14) + 4681; i <= n; i++) {
            bins = addBin(bins, idx++, i);
        }
        return Arrays.copyOf(bins, idx);
    }
    
    private static int[] addBin(int[] bins, int idx, int nbin) {
        if (idx >= bins.length) {
            bins = Arrays.copyOf(bins, bins.length * 2);
        }
        bins[idx] = nbin;
        return bins;
    }
}
