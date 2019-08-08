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

import es.elixir.bsc.ngs.nova.io.DataReaderHelper;
import es.elixir.bsc.ngs.nova.io.DataWriterHelper;
import es.elixir.bsc.ngs.nova.sam.SAMRecord;
import es.elixir.bsc.ngs.nova.sam.SAMTag;
import es.elixir.bsc.ngs.nova.sam.SequenceRecord;
import es.elixir.bsc.ngs.nova.sam.tag.MD;
import es.elixir.bsc.ngs.nova.sam.tag.SAMTagEnum;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.zip.DataFormatException;

/**
 * @author Dmitry Repchevsky
 */

public class BAMRecord extends SAMRecord implements SequenceRecord {
    
    private int refID;
    private int l_seq;
    private int next_refID;
    private short bin;
    
    private byte[] seq;
    private byte[] qual;
    private byte[] auxiliary;
    
    /**
     * Returns BAI index bin as in SAMv1 4.2.1
     * 
     * @return index bin
     */
    public int getBin() {
        if (isUnmappedSegment()) {
            return refID == -1 ? 4680 : Bin.PSEUDO_BIN; // pseudo-bins : reg2bin(-1, 0)
        }

        return BAI.reg2bin(pos - 1, getPositionEnd());
    }

    public int getRefID() {
        return refID;
    }
    
    public void setRefID(final int refID) {
        this.refID = refID;
    }

    public void setNext_refID(final int next_refID) {
        this.next_refID = next_refID;
    }
    
    public int getNextRefID() {
        return next_refID;
    }
    
    @Override
    public String getSequence() {
        return Sequence.parse(seq, l_seq);
    }
    
    @Override
    public void setSequence(final String sequence) {
        l_seq = sequence.length();
        seq = Sequence.parse(sequence);
    }

    private void setSequence(final byte[] sequence, final int l_seq) {
        this.seq = sequence;
        this.l_seq = l_seq;
    }
    
    @Override
    public String getQuality() {
        return Quality.parse(qual);
    }

    @Override
    public void setQuality(final String sequence) {
        this.qual = Quality.parse(sequence);
    }

    private void setQuality(final byte[] quality) {
        this.qual = quality;
    }

    /**
     * <p>
     * Get a standard SAM tag from the record.
     * </p>
     * 
     * Method allows to get a particular standard tag:
     * <br>
     * MD tag = getTagValue(SAMTagEnum.MD)
     * 
     * @param <T> SAM tag type to return (i.e. 'MD')
     * @param tag SAM tag id to return
     * 
     * @return SAM tag or null if no tag found.
     */
    public <T extends SAMTag> T getTag(final SAMTagEnum tag) {
        return (T)getTag(tag.name());
    }
    
    /**
     * <p>
     * Get a SAM tag from the record.
     * </p>
     * 
     * @param tag the name of the SAM tag to get
     * @return the SAM tag or null if no tag found
     */
    public SAMTag getTag(final String tag) {
        final ByteBuffer buf = findTag(tag);
        return buf == null ? null : SAMTagEnum.decode(tag, buf);
    }
    
    /**
     * <p>
     * Set the SAM tag and get back the replaced one.
     * </p>
     * 
     * @param <T> SAM tag type to be set (i.e. 'MD')
     * @param tag SAM tag to be set
     * 
     * @return replaced SAM tag or null in no tag was found
     */
    public <T extends SAMTag> T setTag(final T tag) {
        final SAMTag old_tag = removeTag(tag.getTagName());
        final byte[] data = SAMTagEnum.encode(tag);
        if (auxiliary == null) {
            auxiliary = data;
        } else {
            int currentLength = auxiliary.length;
            auxiliary = Arrays.copyOf(auxiliary, auxiliary.length + data.length);
            System.arraycopy(data, 0, auxiliary, currentLength, data.length);
        }
        return (T)old_tag;
    }
    
    /**
     * <p>
     * Remove a standard SAM tag from the record.
     * </p>
     * 
     * Method allows to remove a particular standard tag:
     * <br>
     * MD tag = removeTag(SAMTagEnum.MD)
     * 
     * @param <T>
     * @param tag SAM tag id to remove
     * 
     * @return SAM tag or null if no tag found.
     */
    public <T extends SAMTag> T removeTag(final SAMTagEnum tag) {
        return (T)removeTag(tag.name());
    }
    
    public SAMTag removeTag(final String tag) {
        final ByteBuffer buf = findTag(tag);
        if (buf != null) {
            final int tag_pos = buf.position() - 2;
            final SAMTag o = SAMTagEnum.decode(tag, buf);
            final byte[] new_auxiliary = new byte[auxiliary.length + tag_pos - buf.position()];
            System.arraycopy(auxiliary, 0, new_auxiliary, 0, tag_pos);
            System.arraycopy(auxiliary, buf.position(), new_auxiliary, tag_pos, auxiliary.length - buf.position());
            auxiliary = new_auxiliary;
            return o;
        }
        return null;
    }
    
    private ByteBuffer findTag(final String tag) {
        if (auxiliary == null) {
            return null;
        }
        
        final ByteBuffer buf = ByteBuffer.wrap(auxiliary).order(ByteOrder.LITTLE_ENDIAN);
        while (buf.hasRemaining()) {
            final char ch0 = (char) buf.get();
            final char ch1 = (char) buf.get();
            if (ch0 == tag.charAt(0) &&
                ch1 == tag.charAt(1)) {
                return buf;
            }
            SAMTagEnum.decode(tag, buf);
        }
        return null;
    }

    /**
     * Reconstructs the reference sequence the read is aligned to.
     * The upper case letters are exact ref matching, while lower case ones 
     * are from the read and may be used to calculate the consent from various 
     * reads.
     * 
     * @return the reconstructed reference sequence
     */
    public String getAlignment() {
        if (cigar == null) {
            return "";
        }

        final StringBuilder sb = new StringBuilder();
        
        final String md = getTag(SAMTagEnum.MD).toString();
        for (int i = 0, idx = 0, n = cigar.length; i < n; i++) {
            
            final long cigar_op = cigar[i];
            final int op = (int)(cigar_op & 0x0F);
            int len = (int)(cigar_op >>> 4);
            switch(op) {
                case 0: // "M" - match or mismatch
                        if (md == null || md.isEmpty()) {
                            // there is no way to know whether it is match or mismatch
                            do {
                                byte ch = idx % 2 == 0 ? Sequence.SEQ[seq[idx++/2] >> 4 & 0x0F] : Sequence.SEQ[seq[idx++/2] & 0x0F];
                                sb.append(Character.toLowerCase(ch));
                            } while (--len > 0);
                            break;
                        }
                case 7: // "=" - seq match
                        do {
                            byte ch = idx % 2 == 0 ? Sequence.SEQ[seq[idx++/2] >> 4 & 0x0F] : Sequence.SEQ[seq[idx++/2] & 0x0F];
                            sb.append((char)ch);
                        } while (--len > 0);
                        break;
                case 5: // "H" - hard clipping
                        break; 
                case 8: // "X" - seq mismatch (substitution)
                        idx += len;
                case 2: // "D" - deletion
                case 3: // "N" - skip
                        do {
                            sb.append('N');
                        } while (--len > 0);
                        break;
                case 1: // "I" - insertion
                case 4: // "S" - soft clipping
                case 6: // "P" - padding
                        idx += len; // skip inserted nucleotides
                        break;
            }
        }
        
        
        if (md != null && md.length() > 0) {
            int idx = 0;
            final Matcher m = MD.PATTERN.matcher(md);
            while (m.find()) {
                final String g1 = m.group(1);
                final String g2 = m.group(2);
                if (g1 != null) {
                    idx += Integer.parseInt(g1);
                }
                if (g2 != null) {
                    final int s = g2.charAt(0) == '^' ? 1 : 0;
                    
                    // MD may go beyond the actual SEQ (i.g. seq is soft clipped)
                    final int n = idx - sb.length() + g2.length();
                    for (int i = s; i < n; i++) {
                        sb.append("N");
                    }

                    for (int i = s; i < g2.length(); i++) {
                        sb.setCharAt(idx++, g2.charAt(i));
                    }                        
                }
            }
        }
        
        return sb.toString();
    }

    public void write(final OutputStream out) throws IOException {
        
        final int l_read_name = qname == null ? 0 : qname.length();
        final int n_cigar_op = cigar == null ? 0 : cigar.length;
        final int l_auxiliary = auxiliary == null ? 0 : auxiliary.length;
        
        final long block_size = 33 + l_read_name + n_cigar_op * 4 + (seq == null ? 0 : seq.length) + l_seq + l_auxiliary;
        
        DataWriterHelper.writeUnsignedInt(out, block_size);
        
        DataWriterHelper.writeUnsignedInt(out, refID);
        DataWriterHelper.writeUnsignedInt(out, pos - 1);
        out.write(l_read_name + 1);
        out.write(mapq);
        DataWriterHelper.writeUnsignedShort(out, getBin());
        DataWriterHelper.writeUnsignedShort(out, n_cigar_op);
        DataWriterHelper.writeUnsignedShort(out, flag);
        DataWriterHelper.writeUnsignedInt(out, l_seq);
        DataWriterHelper.writeUnsignedInt(out, next_refID);
        DataWriterHelper.writeUnsignedInt(out, next_pos - 1);
        DataWriterHelper.writeUnsignedInt(out, tlen);

        if (qname != null) {
            out.write(qname.getBytes(StandardCharsets.US_ASCII));
        }
        out.write(0);
        
        for (int i = 0; i < n_cigar_op; i++) {
            DataWriterHelper.writeUnsignedInt(out, cigar[i]);
        }
        
        if (l_seq > 0) {
            out.write(seq);
            if (qual != null) {
                out.write(qual);
            } else {
                // if there is no qualities defined for the record write '!'
                for (int i = 0; i < l_seq; i++) {
                    out.write(0);
                }
            }
        }
     
        if (l_auxiliary > 0) {
            out.write(auxiliary);
        }
    }

    public static BAMRecord decode(final InputStream in) 
                        throws IOException, DataFormatException {
        
        final BAMRecord record = new BAMRecord();

        final long block_size = DataReaderHelper.readUnsignedInt(in);

        record.refID = (int)DataReaderHelper.readUnsignedInt(in);
        record.pos = (int)DataReaderHelper.readUnsignedInt(in) + 1;
        
        final int l_read_name = in.read() & 0xFF;
        record.mapq = (byte)in.read();
        record.bin = (short)DataReaderHelper.readUnsignedShort(in);
        
        //long flag_nc = DataReader.readUnsignedInt(in);
        final int n_cigar_op = DataReaderHelper.readUnsignedShort(in);
        final int flag = DataReaderHelper.readUnsignedShort(in);
        record.setFlag(flag);
        
        final int l_seq = (int)DataReaderHelper.readUnsignedInt(in);
        
        record.next_refID = (int)DataReaderHelper.readUnsignedInt(in);
        record.next_pos = (int)DataReaderHelper.readUnsignedInt(in) + 1;
        
        record.tlen = (int)DataReaderHelper.readUnsignedInt(in);
        
        StringBuilder read_name = new StringBuilder(l_read_name);
        int ch;
        while ((ch = in.read()) != 0) {
            read_name.append((char)ch);
        }
        record.setQName(read_name.toString());

        
        //int n_cigar_op = (int) (flag_nc & 0xFFFF);
        if (n_cigar_op > 0) {
            final long[] cigar = new long[n_cigar_op];
            for (int i = 0; i < n_cigar_op; i++) {
                final long cigar_op = DataReaderHelper.readUnsignedInt(in);
                cigar[i] = cigar_op;
            }
            record.setCIGAR(cigar);
        }
        
        if (l_seq > 0) {
            final byte[] seq = new byte[(l_seq + 1) / 2];
            for (int i = 0; i < seq.length; i++) {
                seq[i] = (byte)in.read();
            }
            record.setSequence(seq, l_seq);
            
            long sum = 0;
            final byte[] qual = new byte[l_seq];
            for (int i = 0; i < l_seq; i++) {
                qual[i] = (byte)in.read();
                sum += qual[i];
            }
            
            // no need to use memory if all qualities are '!'
            if (sum > 0) {
                record.setQuality(qual);
            }
        }
        
        
        final long read = 33 + read_name.length() + n_cigar_op * 4 + record.seq.length + l_seq;

        int left = (int)(block_size - read);
        if (left > 0) {
            record.auxiliary = new byte[left];
            int idx = 0;
            do {
                idx += in.read(record.auxiliary, idx, left - idx);
            } while (idx < left);
        } else {
            record.auxiliary = null;
        }
        return record;
    }

    public String getGroup() {
        SAMTag tag = getTag(SAMTagEnum.RG);
        return (String) tag.getTagValue();
    }
}
