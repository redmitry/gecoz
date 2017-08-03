package es.elixir.bsc.ngs.nova.fasta;

/**
 *
 * @author Dmitry Repchevsky
 */
public class FastaSequence implements Comparable<FastaSequence> {
    public final String header;
    public final long position;
    public final int length;
    public final byte[] sequence;
    public final boolean multiline;

    public FastaSequence(final String header, final long position, final int length, final boolean multiline) {
        this.header = header;
        this.position = position;
        this.length = length;
        this.sequence = null;
        this.multiline = multiline;
    }
    
    public FastaSequence(final String header, final long position, final byte[] sequence, final boolean multiline) {
        this.header = header;
        this.position = position;
        this.length = sequence.length;
        this.sequence = sequence;
        this.multiline = multiline;
    }

    @Override
    public int compareTo(FastaSequence o) {
        if (this.length != o.length) {
            return this.length > o.length ? -1 : 1;
        }

        return this.header.compareTo(o.header);
    }
}
