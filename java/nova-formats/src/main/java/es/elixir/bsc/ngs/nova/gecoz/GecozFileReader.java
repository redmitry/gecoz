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

import es.elixir.bsc.ngs.nova.algo.ssa.GSSA;
import es.elixir.bsc.ngs.nova.algo.ssa.GSSAIndex;
import es.elixir.bsc.ngs.nova.algo.tree.HSWTShape;
import es.elixir.bsc.ngs.nova.algo.tree.HuffmanShapedWaveletTree;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.StandardOpenOption;
import static java.nio.file.StandardOpenOption.READ;
import java.util.Collections;
import java.util.EnumSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;

/**
 * @author Dmitry Repchevsky
 */

public class GecozFileReader {
    
    private final Map<GecozRefBlockHeader, Long> headers;

    private final FileChannel ref_channel;
    private final FileChannel ssa_channel;
    
    public GecozFileReader(Path path) throws IOException, DataFormatException {
        
        ref_channel = FileChannel.open(path, EnumSet.of(READ));
        
        PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:**.gcz");
        String ssa_fname = path.getFileName().toString();
        if (matcher.matches(path)) {
            ssa_fname = ssa_fname.substring(0, ssa_fname.length() - 3);
        }
        
        Path ssa_path = path.resolveSibling(ssa_fname + "gcx");

        ssa_channel = Files.isReadable(ssa_path) ? FileChannel.open(ssa_path, EnumSet.of(READ)) : null;

        Map<GecozRefBlockHeader, Long> _headers = new LinkedHashMap<>();
        
        long position = 0;
        do {
            InputStream in = Channels.newInputStream(ref_channel);
            GecozRefBlockHeader header = new GecozRefBlockHeader(in);
            _headers.put(header, position);
            position += header.size;
            ref_channel.position(position);
        } while (position < ref_channel.size());
        
        this.headers = Collections.unmodifiableMap(_headers);
    }
    
    public GecozRefBlockHeader findBlockHeader(String header) {
        for (GecozRefBlockHeader hdr : headers.keySet()) {
            if (hdr.findHeader(header) >= 0) {
                return hdr;
            }
        }
        return null;
    }
    
    public Set<GecozRefBlockHeader> getBlockHeaders() {
        return Collections.unmodifiableSet(headers.keySet());
    }

    /**
     * Reads the Succinct Suffix Array from a disk.
     * 
     * @param header the header of the sequence to choose the SSA to read.
     * 
     * @return
     * @throws IOException
     * @throws DataFormatException 
     */
    public GSSA read(GecozRefBlockHeader header) throws IOException, DataFormatException {
        
        Long pos = headers.get(header);
        if (pos == null) {
            return null;
        }
        
        final int hlen = header.getBlockHeaderLength();
        ByteBuffer in = ref_channel.map(FileChannel.MapMode.READ_ONLY, pos + hlen, header.size - hlen);
        in.order(ByteOrder.LITTLE_ENDIAN);
        
        HSWTShape shape = HSWTShape.read(in, header.len);
                
        HuffmanShapedWaveletTree tree = HuffmanShapedWaveletTree.read(shape, in);
        
        if (ssa_channel == null) {
            return new GSSA(tree, null);
        }
        
        // to be able to locate the index in the SSA file we have to recover sampling factor.
        int ssa_headers_length = headers.size() * GecozSSABlockHeader.getBlockHeaderLength();

        // total index data length (minus headers)
        final long ssa_data_length = ssa_channel.size() - ssa_headers_length;

        long test;
        int sampling_factor = -1;
        do {
            test = 0;
            sampling_factor++;
            for (GecozRefBlockHeader h : headers.keySet()) {
                test += GSSAIndex.getIndexSize(h.len, sampling_factor);
            }

        } while (ssa_data_length < test);

        long ssa_pos = 0;
        for (GecozRefBlockHeader h : headers.keySet()) {
            if (header == h) {
                break;
            }
            ssa_pos += GecozSSABlockHeader.getBlockHeaderLength() + GSSAIndex.getIndexSize(h.len, sampling_factor);
        }

        long ssa_size = GSSAIndex.getIndexSize(header.len, sampling_factor);
        ByteBuffer ssa_idx = ssa_channel.map(FileChannel.MapMode.READ_ONLY, ssa_pos, GecozSSABlockHeader.getBlockHeaderLength() + ssa_size);
        ssa_idx.order(ByteOrder.LITTLE_ENDIAN);

        GecozSSABlockHeader ssa_header = new GecozSSABlockHeader(ssa_idx);

        if (header.getHeaderHash() != ssa_header.hash) {
            Logger.getLogger(GecozFileReader.class.getName()).log(Level.SEVERE, "unequal headers");
            throw new DataFormatException("invalid index file");
        }
        if (ssa_header.len != ssa_size) {
            Logger.getLogger(GecozFileReader.class.getName()).log(Level.SEVERE, "unequal header lengths");
            throw new DataFormatException("invalid index file");
        }

        GSSAIndex index = new GSSAIndex(ssa_idx, header.len);

        return new GSSA(tree, index);
    }
    
    @Override
    public String toString() {
        return headers.toString();
    }
    
    public final static boolean checkFormat(Path path) throws IOException {
        try(InputStream in = Files.newInputStream(path, StandardOpenOption.READ)) {
            for (int i = 0, n = GecozRefBlockHeader.MAGIC.length(); i < n; i++) {
                if (GecozRefBlockHeader.MAGIC.charAt(i) != (in.read() & 0xFF)) {
                    return false;
                }
            }
        }
        return true;
    }
}
