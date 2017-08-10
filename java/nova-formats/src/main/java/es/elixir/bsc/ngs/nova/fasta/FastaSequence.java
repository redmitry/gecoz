package es.elixir.bsc.ngs.nova.fasta;

/**
 *
 * @author Dmitry Repchevsky
 */
public class FastaSequence extends TFastaSequence {
    public final long position;
    public final byte[] sequence;

    public FastaSequence(final String header, final long position, final int length, final boolean multiline) {
        super(header, length, multiline);
        this.position = position;
        this.sequence = null;
    }
    
    public FastaSequence(final String header, final long position, final byte[] sequence, final boolean multiline) {
        super(header, sequence.length, multiline);
        this.position = position;
        this.sequence = sequence;
    }
}
