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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.atomic.AtomicLong;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;

/**
 * @author Dmitry Repchevsky
 */

public class GecoMatch {

    /**
     * Searches the pattern in the GecoZ file
     * 
     * @param ipath the path for the GecoZ file
     * @param header the sequence header where to search the pattern (or null)
     * @param pattern the pattern to search
     */
    static void match(Path ipath, String header, String pattern) {
        match(ipath, header, pattern, true);
    }
    
    static void count(Path ipath, String header, String pattern) {
        match(ipath, header, pattern, false);
    }
    
    private static void match(Path ipath, String header, String pattern, boolean match) {
        try {
            if (!Files.exists(ipath) || Files.isDirectory(ipath)) {
                Logger.getLogger(GecoMatch.class.getName()).log(Level.SEVERE, "no gecoz file found: {0}\n", ipath);
                System.exit(1); 
            }

            if (!GecozFileReader.checkFormat(ipath)) {
                Logger.getLogger(GecoMatch.class.getName()).log(Level.SEVERE, "invalid gecoz file format: {0}\n", ipath);
                System.exit(1); 
            }

            GecozFileReader reader = new GecozFileReader(ipath);
            
            if (header != null) {
                GecozRefBlockHeader bheader = reader.findBlockHeader(header);
                if (bheader == null) {
                    Logger.getLogger(GecoMatch.class.getName()).log(Level.SEVERE, "no sequence found: {0}\n", header);
                    System.exit(1);
                }
                
                GSSA ssa = reader.read(bheader);
                if (ssa == null) {
                    Logger.getLogger(GecoMatch.class.getName()).log(Level.SEVERE, "no sequence found: {0}\n", header);
                    System.exit(1);
                }

                final int nstr = bheader.findHeader(header);
                
                final long t1 = System.nanoTime();
                long[][] res = ssa.find(pattern.getBytes("UTF8"));
                final long t2 = System.nanoTime();

                if (res != null && res.length > 0 && res[nstr] != null && res[nstr].length > 0) {
                    System.out.println(">" + header + " found : " + res[nstr].length);
                    if (match) {
                        for (int i = 0; i < res[nstr].length; i++) {
                            System.out.println(res[nstr][i]);
                        }
                    }
                }
                Logger.getLogger(GecoMatch.class.getName()).log(Level.INFO, "finished in {0} ms.", (t2 - t1)/1000000);
            } else {
                long count = 0; // total matches counter
                
                final AtomicLong time = new AtomicLong();  // total time spent
                
                for (GecozRefBlockHeader bheader : reader.getBlockHeaders()) {
                    final long t = System.nanoTime();
                    
                    GSSA ssa = reader.read(bheader);
                    if (ssa == null) {
                        Logger.getLogger(GecoMatch.class.getName()).log(Level.SEVERE, "no block found: {0} skipping...\n", bheader.len);
                        continue;
                    }

                    long[][] res = ssa.find(pattern.getBytes("UTF8"));

                    time.addAndGet(System.nanoTime() - t);

                    if (res != null && res.length > 0) {
                        // print results
                        count += print(bheader.headers, res, match);
                    }
                }
                
                Logger.getLogger(GecoMatch.class.getName()).log(Level.INFO, "total found: {0}\n", count);
                Logger.getLogger(GecoMatch.class.getName()).log(Level.INFO, "finished in {0} ms.\n", time.get()/1000000);
            }
        } catch(IOException | DataFormatException ex) {
            Logger.getLogger(GecoMatch.class.getName()).log(Level.SEVERE, "error reading a file: {0}\n", ipath);
            Logger.getLogger(GecoMatch.class.getName()).log(Level.SEVERE, ex.getMessage());
            System.exit(1);
        }
    }
    
    private static long print(String[] headers, long[][] res, boolean match) {
        long count = 0;
        for (int i = 0, n = res.length; i < n; i++) {
            if (res[i] != null && res[i].length > 0) {
                count += res[i].length;
                System.out.println(">" + headers[i] + " found : " + res[i].length);
                if (match) {
                    for (int j = 0; j < res[i].length; j++) {
                        System.out.println(res[i][j]);
                    }
                }
            }
        }
        return count;
    }
}
