package es.elixir.bsc.ngs.nova.fasta;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Dmitry Repchevsky
 */

public class FastaIterator implements Iterator<FastaSequence>, AutoCloseable {
    
    private int ch;
    private long position;
    private String header;
    private final InputStream in;
    private final ByteArrayOutputStream out;
    
    public FastaIterator(InputStream in, boolean lazy) {
        this.in = in;
        out = lazy ? null : new ByteArrayOutputStream();
    }
    
    protected FastaIterator(InputStream in, FastaSequence seq, boolean lazy) throws IOException {
        ch = '\r';
        position = seq.position;
        header = seq.header;
        out = lazy ? null : new ByteArrayOutputStream();
        
        this.in = in;
        long toskip = position;
        while ((toskip -= in.skip(toskip)) > 0); // put stream to the sequence position
    }
    
    @Override
    public boolean hasNext() {
        if (header != null) {
            return true;
        }

        try {
            while (ch >= 0 && ch != '>' && ch != '@') {
                ch = in.read();
                position++;
            }

            if (ch < 0) {
                return false;
            }

            StringBuilder sb = new StringBuilder();
            while ((ch = in.read()) >= 0 && ch != '\n') {
                position++;
                if (ch != '\r') {
                    sb.append((char)ch);
                }
            }
            position++;
            header = sb.toString();
        } catch (IOException ex) {
            return false;
        }

        return true; // return !header.isEmpty();
    }

    @Override
    public FastaSequence next() {
        if (!hasNext()) {
            return null;
        }

        if (out != null) {
            out.reset();
        }

        int lines = 0;
        int length = 0;
        long posnew = position;
        try {
            do {
                if (ch >= 0 && ch != '\r' && ch != '\n') {
                    lines++;
                    do {
                        if (out != null) {
                            out.write(ch);
                        }
                        posnew++;
                        length++;
                    } while ((ch = in.read()) >= 0 && ch != '\r' && ch != '\n');
                }
                posnew++;
            } while((ch = in.read()) >= 0 && ch != '>' && ch !='@' && ch !='+');
            
            if (ch == '+') {
                // skip qualities
                int qlines = -1;
                int qlength = 0;
                do {
                    while ((ch = in.read()) >= 0 && ch != '\r' && ch != '\n') {
                        qlength++;
                        posnew++;
                    }
                    posnew++;
                    qlines++;
                } while(qlength < length && qlines < lines);
                
                if (length != qlength) {
                    Logger.getLogger(FastaIterator.class.getSimpleName()).log(Level.WARNING, "different sequence and qualities lengths ''{0}'' ({1} bytes)\n", new Object[]{length, qlength});
                }
            }
        } catch (IOException ex) {
            return null;
        }

        final FastaSequence sequence = out != null ? 
                new FastaSequence(header, position, out.toByteArray(), lines > 1) :
                new FastaSequence(header, position, length, lines > 1);
        
        header = null;
        position = posnew;
        return sequence;
    }

    public String getHeader() {
        return header;
    }
    
    @Override
    public void close() throws Exception {
        in.close();
    }
}
