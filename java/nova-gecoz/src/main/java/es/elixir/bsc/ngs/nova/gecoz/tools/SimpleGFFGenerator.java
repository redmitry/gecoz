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

import es.elixir.bsc.ngs.nova.algo.ssa.GSSA;
import es.elixir.bsc.ngs.nova.gecoz.GecozFileReader;
import es.elixir.bsc.ngs.nova.gecoz.GecozRefBlockHeader;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Set;
import java.util.zip.DataFormatException;

/**
 * @author Dmitry Repchevsky
 */

public class SimpleGFFGenerator {
    
    public static void search(Path ref, Path fasta) {

        try {            
            GecozFileReader reader = new GecozFileReader(ref);
            
            Set<GecozRefBlockHeader> set = reader.getBlockHeaders();
            GecozRefBlockHeader[] bheaders = new GecozRefBlockHeader[set.size()];
            bheaders = set.toArray(bheaders);
            GSSA[] gssa = new GSSA[bheaders.length];
            
            for (int i = 0, n = gssa.length; i < n; i++) {
                gssa[i] = reader.read(bheaders[i]);
            }
            
            ByteArrayOutputStream o = new ByteArrayOutputStream();
            try (BufferedReader fr = Files.newBufferedReader(fasta)) {
                String line;
                String header = null;
                while ((line = fr.readLine()) != null) {
                    if (line.startsWith(">") || line.startsWith("@")) {
                        if (header != null && o.size() != 0) {
                            final byte[] arr = o.toByteArray();
                            search(gssa, bheaders, header, arr);
                        }
                        header = line.substring(1);
                        o.reset();
                    } else if (line.startsWith("+")) {
                        if (header != null && o.size() != 0) {
                            final byte[] arr = o.toByteArray();
                            search(gssa, bheaders, header, arr);
                        }
                        header = null;
                        o.reset();
                    } else if (header != null) {
                        o.write(line.getBytes("UTF8"));
                    }
                }

                if (header != null && o.size() != 0) {
                    final byte[] arr = o.toByteArray();
                    search(gssa, bheaders, header, arr);
                }
            }
        } catch(IOException | DataFormatException ex) {
            System.err.println("error reading file: " + ref);
            ex.printStackTrace(System.err);
            System.exit(1);
        }
    }

    private static void search(GSSA[] gssa, GecozRefBlockHeader[] bheaders, String header, byte[] seq) throws IOException, DataFormatException {
        
        for (int i = 0, n = seq.length; i < n; i++) {
            if (seq[i] == 'U') {
                seq[i] = 'T';
            }
        }

        search(gssa, bheaders, header, seq, false);
        
        for (int i = 0, n = seq.length - 1; i <= n;) {
            final byte b = reverse(seq[i]);
            seq[i++] = reverse(seq[n]);
            seq[n--] = b;
        }

        search(gssa, bheaders, header, seq, true);
    }

    private static byte reverse(byte b) {
        switch (b) {
            case 'A' : return 'T';
            case 'T' : return 'A';
            case 'C' : return 'G';
            case 'G' : return 'C';
            default  : return b;
        }
    }

    private static void search(GSSA[] gssa, GecozRefBlockHeader[] bheaders, String header, byte[] seq, boolean reverse) throws IOException, DataFormatException {

        for (int i = 0; i < gssa.length; i++) {
            
            if (gssa[i] != null) {
                long[][] res = gssa[i].find(seq);

                if (res != null) {
                    for (int j = 0; j < res.length; j++) {
                        if (res[j] != null) {
                            for (int k = 0; k < res[j].length; k++) {
                                System.out.print(bheaders[i].headers[j]);
                                System.out.print('\t');
                                System.out.print("gecotools");
                                System.out.print('\t');
                                System.out.print("dna");
                                System.out.print('\t');
                                System.out.print(res[j][k] + 1);
                                System.out.print('\t');
                                System.out.print(res[j][k] + seq.length);
                                System.out.print("\t1.000\t");
                                System.out.print(reverse ? '-' : '+');
                                System.out.print("\t.\t");

                                String[] h = header.split("\\|");
                                if (h.length > 0) {
                                    System.out.print("ID=");
                                    System.out.print(h[0]);
                                }
                                for(int l = 1; l < h.length; l++) {
                                    System.out.print(";Note=");
                                    System.out.print(h[l]);
                                }
                                System.out.print('\n');
                            }
                        }
                    }
                }
            }
        }
    }
}
