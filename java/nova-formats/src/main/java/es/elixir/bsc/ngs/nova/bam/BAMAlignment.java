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
import es.elixir.bsc.ngs.nova.sam.MD;
import es.elixir.bsc.ngs.nova.sam.SAMRecord;
import es.elixir.bsc.ngs.nova.sam.SAMTag;
import es.elixir.bsc.ngs.nova.sam.SequenceRecord;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.regex.Matcher;
import java.util.zip.DataFormatException;

/**
 * @author Dmitry Repchevsky
 */

public class BAMAlignment extends SAMRecord implements SequenceRecord {
    
    private final char[] VAL_TYPE = {'A','c','C','s','S','i','I','f','Z','H','B'};
    
    private int refID;
    private int l_seq;
    private int next_refID;
    private int next_pos;
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
        if (pos == 0) {
            return 4680; // reg2bin(-1, 0)
        }
        return BAI.reg2bin(pos - 1, getPositionEnd());
    }

    public int getRefID() {
        return refID;
    }
    
    public int getNextRefID() {
        return next_refID;
    }
    
    public int getNextPosD() {
        return next_pos;
    }
    
    @Override
    public String getSequence() {
        return Sequence.parse(seq, l_seq);
    }
    
    @Override
    public void setSequence(String sequence) {
        l_seq = sequence.length();
        seq = Sequence.parse(sequence);
    }

    private void setSequence(byte[] sequence, int l_seq) {
        this.seq = sequence;
        this.l_seq = l_seq;
    }
    
    @Override
    public void setQuality(String sequence) {
        this.qual = Quality.parse(sequence);
    }

    private void setQuality(byte[] quality) {
        this.qual = quality;
    }

    @Override
    public String getQuality() {
        return Quality.parse(qual);
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

        StringBuilder sb = new StringBuilder();
        
        String md = (String)getTagValue(SAMTag.MD);
        for (int i = 0, idx = 0, n = cigar.length; i < n; i++) {
            
            final long cigar_op = cigar[i];
            final int op = (int)(cigar_op & 0x0F);
            int len = (int)(cigar_op >>> 4);
            switch(op) {
                case 0: // "M" - match or mismatch
                        if (md == null || md.isEmpty()) {
                            // there is no way to know whether it is match or mismatch
                            do {
                                char ch = idx % 2 == 0 ? Sequence.SEQ[seq[idx++/2] >> 4 & 0x0F] : Sequence.SEQ[seq[idx++/2] & 0x0F];
                                sb.append(Character.toLowerCase(ch));
                            } while (--len > 0);
                            break;
                        }
                case 7: // "=" - seq match
                        do {
                            char ch = idx % 2 == 0 ? Sequence.SEQ[seq[idx++/2] >> 4 & 0x0F] : Sequence.SEQ[seq[idx++/2] & 0x0F];
                            sb.append(ch);
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
            Matcher m = MD.PATTERN.matcher(md);
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
    
    public Object getTagValue(final SAMTag tag) {
        if (auxiliary == null) {
            return null;
        }
        ByteBuffer buf = ByteBuffer.wrap(auxiliary).order(ByteOrder.LITTLE_ENDIAN);
        
        while (buf.hasRemaining()) {
            //String nm = new String(auxiliary, ++idx, 2);
            char ch0 = (char) buf.get();
            char ch1 = (char) buf.get();
            
            Object o = SAMTag.decode(buf);
            if (ch0 == tag.name().charAt(0) &&
                ch1 == tag.name().charAt(1)) {
                return o;
            }
        }
        
        return null;
    }
    
    public void write() {
        
    }

    public static BAMAlignment decode(final InputStream in) 
                        throws IOException, DataFormatException {
        
        BAMAlignment alignment = new BAMAlignment();
        final long block_size = DataReaderHelper.readUnsignedInt(in);
        alignment.refID = (int)DataReaderHelper.readUnsignedInt(in);
        
        final int pos = (int)DataReaderHelper.readUnsignedInt(in);
        alignment.setPositionStart(pos + 1);
        
        final int l_read_name = in.read() & 0xFF;
        alignment.mapq = (byte)in.read();
        alignment.bin = (short)DataReaderHelper.readUnsignedShort(in);
        
        //long flag_nc = DataReader.readUnsignedInt(in);
        int n_cigar_op = DataReaderHelper.readUnsignedShort(in);
        final int flag = DataReaderHelper.readUnsignedShort(in);
        alignment.setFlag(flag);
        
        final int l_seq = (int)DataReaderHelper.readUnsignedInt(in);
        
        alignment.next_refID = (int)DataReaderHelper.readUnsignedInt(in);
        alignment.next_pos = (int)DataReaderHelper.readUnsignedInt(in);
        
        alignment.tlen = (int)DataReaderHelper.readUnsignedInt(in);
        
        StringBuilder read_name = new StringBuilder(l_read_name);
        int ch;
        while ((ch = in.read()) != 0) {
            read_name.append((char)ch);
        }
        alignment.setQName(read_name.toString());

        
        //int n_cigar_op = (int) (flag_nc & 0xFFFF);
        if (n_cigar_op > 0) {
            final long[] cigar = new long[n_cigar_op];
            for (int i = 0; i < n_cigar_op; i++) {
                final long cigar_op = DataReaderHelper.readUnsignedInt(in);
                cigar[i] = cigar_op;
            }
            alignment.setCIGAR(cigar);
        }
        
        if (l_seq > 0) {
            final byte[] seq = new byte[(l_seq + 1) / 2];
            for (int i = 0; i < seq.length; i++) {
                seq[i] = (byte)in.read();
            }
            alignment.setSequence(seq, l_seq);
            
            final byte[] qual = new byte[l_seq];
            for (int i = 0; i < l_seq; i++) {
                qual[i] = (byte)in.read();
            }
            alignment.setQuality(qual);
        }
        
        
        final long read = 32 + read_name.length() + 1 + n_cigar_op * 4 + alignment.seq.length + l_seq;

        int left = (int)(block_size - read);
        if (left > 0) {
            alignment.auxiliary = new byte[left];
            int idx = 0;
            do {
                idx += in.read(alignment.auxiliary, idx, left - idx);
            } while (idx < left);
        } else {
            alignment.auxiliary = null;
        }
        return alignment;
    }
}
