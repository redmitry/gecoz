package es.elixir.bsc.ngs.nova.fasta;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

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
