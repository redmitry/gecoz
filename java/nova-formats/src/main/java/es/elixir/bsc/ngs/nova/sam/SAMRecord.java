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

package es.elixir.bsc.ngs.nova.sam;

/**
 * @author Dmitry Repchevsky
 */

public abstract class SAMRecord implements SequenceRecord {
    public final static int UNMAPPED_SEGMENT = 0x04;
    
    protected int pos;
    protected int flag;
    protected byte mapq;
    protected int tlen;
    protected String qname;
    protected String rname;
    protected int end_pos;
    protected long[] cigar;
    
    public byte getMappingQuality() {
        return mapq;
    }

    public void setMappingQuality(final byte mapq) {
        this.mapq = mapq;
    }
    
    public int getTemplateLength() {
        return tlen;
    }
    
    public void setTemplateLength(final int tlen) {
        this.tlen = tlen;
    }

    @Override
    public String getQName() {
        return qname;
    }
    
    @Override
    public void setQName(final String qname) {
        this.qname = qname;
    }
    
    public String getRName() {
        return rname;
    }
    
    public void setRName(final String rname) {
        this.rname = rname;
    }

    @Override
    public String getCIGAR() {
        return cigar != null ? CIGAR.parse(cigar) : "*"; // set '*' if unavailable
    }
    
    @Override
    public void setCIGAR(final String cigar) {
        this.cigar = CIGAR.encode(cigar);
    }

    protected void setCIGAR(final long[] cigar) {
        this.cigar = cigar;
    }

    @Override
    public int getPositionStart() {
        return pos;
    }
    
    @Override
    public void setPositionStart(final int pos) {
        this.pos = pos;
    }
    
    @Override
    public int getPositionEnd() {
        if(isUnmappedSegment()) {
            // For unmapped reads treat the alignment as being length one.
            return 0;
        }
        
        if (end_pos == 0 && cigar != null) {
            end_pos = pos + CIGAR.getLength(cigar) - 1;
        }
        
        return end_pos;
    }

    protected void setPositionEnd(final int end_pos) {
        this.end_pos = end_pos;
    }

    @Override
    public int getFlag() {
        return flag;
    }
    
    @Override
    public void setFlag(final int flag) {
        this.flag = flag;
    }
    
    public boolean isUnmappedSegment() {
        return (flag & UNMAPPED_SEGMENT) != 0;
    }
}
