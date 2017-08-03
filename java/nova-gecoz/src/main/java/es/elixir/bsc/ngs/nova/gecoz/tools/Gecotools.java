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

import es.elixir.bsc.ngs.nova.gecoz.GecozFileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.MessageFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

/**
 * @author Dmitry Repchevsky
 */

public class Gecotools {

    private final static String HELP = 
            "gecotools -i file [optional params]\n\n" +
            "parameters:\n\n" +
            "-h (--help)           - this help message\n" +
            "-i (--input)          - either *.fa or *.gcz\n" +
            "-o [header][from][to] - depends on the input parameters\n" +
            "                        (*.fa -> *.gcz, *.gcz -> *.fa, *gcz -> *.seq)\n" +
            "-c [header] 'string'  - count string occurrences in the *.gcz file\n" +
            "-s [header] 'string'  - search string in the *.gcz file\n" +
            "-t                    - use n threads \n" +
            "-v [level]            - verbose (default = WARNING) \n\n" +
            "examples:\n\n" +
            ">java -jar gecotools.jar -t 8 -i hg38.fa -o hg38.gcz\n" +
            ">java -jar gecotools.jar -i hg38.gcz -o hg38.fasta\n" +
            ">java -jar gecotools.jar -i hg38.gcz -o chr15.seq chr15\n" +
            ">java -jar gecotools.jar -i hg38.gcz -c ATTAACCCATGAAAA\n" +
            ">java -jar gecotools.jar -i hg38.gcz -s ATTAACCCATGAAAA\n" +
            ">java -jar gecotools.jar -i hg38.gcz -s chr11 ATTAACCCATGAAAA\n";

    /**
     * 
     * @param args 
     */
    public static void main(String[] args) {
        Map<String, List<String>> params = parameters(args);
        
        if (params.isEmpty() || 
            params.get("-h") != null ||
            params.get("--help") != null) {
            System.out.println(HELP);
            System.exit(0);            
        }
        
        setVerbosity(params.get("-v"));
        
        List<String> in = params.get("-i");
        if (in == null) {
            in = params.get("--input");
        } else if (params.containsKey("--input")) {
            System.err.println("only one of the forms may be used - either '-i' or '--input'");
            System.exit(1);
        }
        if (in == null) {
            System.err.println("no input file specified");
            System.exit(1);
        } else if (in.size() > 1) {
            System.err.println("more than one input files specified");
            System.exit(1);            
        }

        Path ipath = Paths.get(in.get(0));

        if (params.get("-o") != null) {
            out(ipath, params);
        } else if (params.get("-s") != null) {
            search(ipath, params);
        } else if (params.get("-c") != null) {
            count(ipath, params);
        }
    }
    
    private static void out(Path ipath, Map<String, List<String>> params) {
        List<String> out = params.get("-o");
        if (out.isEmpty()) {
            System.err.println("no output file specified.");
            System.exit(1);
        }
        
        List<String> threads = params.get("-t");
        
        final int th = threads == null || threads.isEmpty() ? 1 : Integer.valueOf(threads.get(0));
        
        Path opath = Paths.get(out.get(0));
        
        try {
            if (GecozFileReader.checkFormat(ipath)) {
                if (out.size() > 1) {
                    // extract sequence by header and positions
                    final String header = out.get(1);
                    final int from = out.size() > 2 ? Integer.valueOf(out.get(2)) : 0;
                    final int to = out.size() > 3 ? Integer.valueOf(out.get(3)) : Integer.MAX_VALUE;
                    
                    GecoRead.sequence(ipath, header, from, to, opath);
                } else {
                    // recover entire fasta file
                    GecoRead.fasta(ipath, opath);
                }
            } else {
                // compress and (optionally) generate index file
                List<String> idx = params.get("-idx");
                Path xpath = idx == null || idx.isEmpty() ? null : Paths.get(idx.get(0));
                
                GecoIndex.index(ipath, opath, xpath, 32, th);
            }
        } catch(IOException ex) {
            System.err.println("error reading a file: " + ipath);
            System.exit(1);
        }
    }
    
    private static void search(Path ref, Map<String, List<String>> params) {
        List<String> search = params.get("-s");
        if (search.size() < 1) {
            System.err.println("no search string/filename specified.");
            System.exit(1);   
        }
        
        final String pattern = search.size() == 1 ? search.get(0): search.get(1);

        if (search.size() == 1) {
            Path fasta = Paths.get(pattern);
            if (Files.isRegularFile(fasta)) {
                SimpleGFFGenerator.search(ref, fasta);
                return;
            }
        }
        
        GecoMatch.match(ref, search.size() > 1 ? search.get(0) : null, pattern);

    }

    private static void count(Path ipath, Map<String, List<String>> params) {
        List<String> count = params.get("-c");
        
        if (count.size() < 1) {
            System.err.println("no search string specified.");
            System.exit(1);   
        }
        
        GecoMatch.count(ipath, 
                        count.size() > 1 ? count.get(0) : null, 
                        count.size() == 1 ? count.get(0): count.get(1));
    }
    
    private static void setVerbosity(List<String> verbosity) {
        Level level;
        try {
            level = verbosity == null || verbosity.isEmpty() ? Level.WARNING : 
                    Level.parse(verbosity.get(0));
        } catch(IllegalArgumentException ex) {
            level = Level.WARNING;
        }

        Formatter formatter = new Formatter() {
            @Override
            public String format(LogRecord record) {
                return MessageFormat.format(record.getMessage(), record.getParameters());
            }
        };

        Logger log = LogManager.getLogManager().getLogger("");
        log.setLevel(level);
        for (Handler h : log.getHandlers()) {
            h.setLevel(level);
            h.setFormatter(formatter);
        }
    }

    private static Map<String, List<String>> parameters(String[] args) {
        TreeMap<String, List<String>> parameters = new TreeMap();        
        List<String> values = null;
        for (String arg : args) {
            switch(arg) {
                case "-h":
                case "--help":
                case "-i":
                case "--input":
                case "-idx":
                case "--index":
                case "-s":
                case "--search":
                case "-c":
                case "--count":
                case "-a":
                case "--align":
                case "-t":
                case "--threads":
                case "-v":
                case "--verbose":
                case "-o":
                case "--output": values = parameters.get(arg);
                              if (values == null) {
                                  values = new ArrayList(); 
                                  parameters.put(arg, values);
                              }
                              break;
                default: if (values != null) {
                    values.add(arg);
                }
            }
        }
        return parameters;
    }
}
