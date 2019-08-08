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

package es.elixir.bsc.ngs.nova.sam.header;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Dmitry Repchevsky
 */

public class SAMHeader {
    protected HeaderLine header_line;
    protected List<ReferenceLine> references;
    protected List<ReadGroupLine> readGroups;
    protected List<ProgramLine> programs;
    protected List<CommentLine> comments;

    protected SAMHeader() {}

    public SAMHeader(final HeaderLine header) {
        this.header_line = header;
    }

    public SAMHeader(final SAMHeader header) {
        this.header_line = header.header_line;
        this.references = header.references;
        this.readGroups = header.readGroups;
        this.programs = header.programs;
        this.comments = header.comments;        
    }
    
    public SAMHeader(
            final HeaderLine header,
            final List<ReferenceLine> references,
            final List<ReadGroupLine> readGroups,
            final List<ProgramLine> programs,
            final List<CommentLine> comments) {

        this.header_line = header;
        this.references = references;
        this.readGroups = readGroups;
        this.programs = programs;
        this.comments = comments;
    }
    
    public SAMHeader(final String header) {
        final String[] lines = header.split("\n");
        if(!lines[0].startsWith(HeaderLine.TAG)) {
            throw new IllegalArgumentException();
        }
        header_line = new HeaderLine(lines[0].substring(4));

        for(int i = 1; i < lines.length; i++) {
            switch (lines[i].substring(0, 3)) {
                case ReferenceLine.TAG:
                    getReferences().add(new ReferenceLine(lines[i].substring(4)));
                    break;
                case ReadGroupLine.TAG:
                    getReadGroups().add(new ReadGroupLine(lines[i].substring(4)));
                    break;
                case ProgramLine.TAG:
                    getPrograms().add(new ProgramLine(lines[i].substring(4)));
                    break;
                case CommentLine.TAG:
                    getComments().add(new CommentLine(lines[i].substring(4)));
                    break;
                default:
                    throw new IllegalArgumentException();
            }
        }
    }

    public HeaderLine getHeaderLine() {
        return header_line;
    }

    public void setHeaderLine(final HeaderLine header_line) {
        this.header_line = header_line;
    }

    public final List<ReferenceLine> getReferences() {
        if (references == null) {
            references = new ArrayList<>();
        }
        return references;
    }

    public final List<ReadGroupLine> getReadGroups() {
        if (readGroups == null) {
            readGroups = new ArrayList<>();
        }
        return readGroups;
    }

    public final List<ProgramLine> getPrograms() {
        if (programs == null) {
            programs = new ArrayList<>();
        }        
        return programs;
    }

    public final List<CommentLine> getComments() {
        if (comments == null) {
            comments = new ArrayList<>();
        }   
        return comments;
    }
    
    public final String[] getGroupIds() {
        String[] groupIds = new String[getReadGroups().size()];
        for(int i = 0; i < groupIds.length; i++){
            groupIds[i] = readGroups.get(i).readGroupId;
        }
        return groupIds;
    }

    public void write(final PrintStream out) throws IOException {

        if (header_line != null) {
            header_line.write(out);
            out.append('\n');
        }

        if (references != null) {
            for(ReferenceLine line : references) {
                line.write(out);
                out.append('\n');
            }
        }
        if (readGroups != null) {
            for(ReadGroupLine line : readGroups) {
                line.write(out);
                out.append('\n');
            }
        }
        if (programs != null) {
            for(ProgramLine line : programs) {
                line.write(out);
                out.append('\n');
            }
        }
        if (comments != null) {
            for(CommentLine line : comments) {
                line.write(out);
                out.append('\n');
            }
        }
    }

    @Override
    public String toString() {
        try (ByteArrayOutputStream out = new ByteArrayOutputStream()) {
            write(new PrintStream(out, false, StandardCharsets.US_ASCII.name()));
            return out.toString(StandardCharsets.US_ASCII.name());
        } catch (IOException ex) {}
        return "";
    }    
}
