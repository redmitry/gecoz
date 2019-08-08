/**
 * *****************************************************************************
 * Copyright (C) 2019 Spanish National Bioinformatics Institute (INB) and
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
    public final static int HAS_MULTIPLE_SEGMENTS = 0x01;
    public final static int SEGMENT_PROPERLY_ALIGNED = 0x02;
    public final static int UNMAPPED_SEGMENT = 0x04;
    public final static int NEXT_SEGMENT_UNMAPPED = 0x08;
    public final static int REVERSE_COMPLEMENTED = 0x10;
    public final static int NEXT_SEGMENT_REVERSE_COMPLEMENTED = 0x20;
    public final static int FIRST_SEGMENT = 0x40;
    public final static int LAST_SEGMENT = 0x80;
    public final static int SECONDARY_ALIGNMENT = 0x100;
    public final static int NOT_PASSING_QUALITY = 0x200;
    public final static int OPTICAL_DUPLICATE = 0x400;
    public final static int SUPPLEMENTARY_ALIGNMENT = 0x800;

    protected int pos;
    protected int next_pos;
    protected int flag;
    protected byte mapq;
    protected int tlen;
    protected String qname;
    protected String rname;
    protected String rname_next;
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

    public String getRNameNext() {
        return rname_next;
    }

    public void setRNameNext(final String rname_next) {
        this.rname_next = rname_next;
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

    public int getNextPositionStart() {
        return next_pos;
    }

    public void setNextPositionStart(final int next_pos) {
        this.next_pos = next_pos;
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

    private void setFlag(final int mask, final boolean value) {
        flag = value ? flag |= mask : flag & ~mask;
    }

    public boolean hasMultipleSegments() {
        return (flag & HAS_MULTIPLE_SEGMENTS) != 0;
    }

    public void setHasMultipleSegments(final boolean value) {
        setFlag(HAS_MULTIPLE_SEGMENTS, value);
    }

    public boolean isSegmentProperlyAligned() {
        return (flag & SEGMENT_PROPERLY_ALIGNED) != 0;
    }

    public void setSegmentProperlyAligned(final boolean value) {
        setFlag(SEGMENT_PROPERLY_ALIGNED, value);
    }

    public boolean isUnmappedSegment() {
        return (flag & UNMAPPED_SEGMENT) != 0;
    }

    public void setUnmappedSegment(final boolean value) {
        setFlag(UNMAPPED_SEGMENT, value);
    }

    public boolean isNextSegmentUnmapped() {
        return (flag & NEXT_SEGMENT_UNMAPPED) != 0;
    }

    public void setNextSegmentUnmapped(final boolean value) {
        setFlag(NEXT_SEGMENT_UNMAPPED, value);
    }

    public boolean isReverseComplemented() {
        return (flag & REVERSE_COMPLEMENTED) != 0;
    }

    public void setReverseComplemented(final boolean value) {
        setFlag(REVERSE_COMPLEMENTED, value);
    }

    public boolean isNextSegmentReverseComplemented() {
        return (flag & NEXT_SEGMENT_REVERSE_COMPLEMENTED) != 0;
    }

    public void setNextSegmentReverseComplemented(final boolean value) {
        setFlag(NEXT_SEGMENT_REVERSE_COMPLEMENTED, value);
    }

    public boolean isFirstSegment() {
        return (flag & FIRST_SEGMENT) != 0;
    }

    public void setFirstSegment(final boolean value) {
        setFlag(FIRST_SEGMENT, value);
    }

    public boolean isLastSegment() {
        return (flag & LAST_SEGMENT) != 0;
    }

    public void setLastSegment(final boolean value) {
        setFlag(LAST_SEGMENT, value);
    }

    public boolean isSecondaryAlignment() {
        return (flag & SECONDARY_ALIGNMENT) != 0;
    }

    public void setSecondaryAlignment(final boolean value) {
        setFlag(SECONDARY_ALIGNMENT, value);
    }

    public boolean isNotPassingQuality() {
        return (flag & NOT_PASSING_QUALITY) != 0;
    }

    public void setNotPassingQuality(final boolean value) {
        setFlag(NOT_PASSING_QUALITY, value);
    }

    public boolean isOpticalDuplicate() {
        return (flag & OPTICAL_DUPLICATE) != 0;
    }

    public void setOpticalDuplicate(final boolean value) {
        setFlag(OPTICAL_DUPLICATE, value);
    }

    public boolean isSupplementaryAlignment() {
        return (flag & SUPPLEMENTARY_ALIGNMENT) != 0;
    }

    public void setSupplementaryAlignment(final boolean value) {
        setFlag(SUPPLEMENTARY_ALIGNMENT, value);
    }
}
