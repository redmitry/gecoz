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

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * @author Dmitry Repchevsky
 */

public class ProgramLine extends AbstractHeaderLine{
    
    public final static String TAG = "@PG";
    
    public final static String[] TAGS = {"ID", "PN", "CL", "PP", "DS", "VN"};
    
    public final String id;
    public final String programName;
    public final String commandLine;
    public final String previousProgramId;
    public final String description;
    public final String programVersion;

    public ProgramLine(
            final String id,
            final String programName,
            final String commandLine,
            final String previousProgramId,
            final String description,
            final String programVersion) {
        this.id = id;
        this.programName = programName;
        this.commandLine = commandLine;
        this.previousProgramId = previousProgramId;
        this.description = description;
        this.programVersion = programVersion;
    }
            
    public ProgramLine(final String line) {
        
        final String[] tags = parseHeaderLine(line, Arrays.copyOf(TAGS, TAGS.length));
        
        id = tags[0];
        programName = tags[1];
        commandLine = tags[2];
        previousProgramId = tags[3];
        description = tags[4];
        programVersion = tags[5];
    }

    @Override
    public void write(final PrintStream out) throws IOException {
        out.append(TAG);
        
        if (id != null) {
            out.append("\tID:").append(id);
        }
        if(programName != null){
            out.append("\tPN:").append(programName);
        }
        if(previousProgramId != null){
            out.append("\tPP:").append(previousProgramId);
        }
        if(description != null){
            out.append("\tDS:").append(description);
        }
        if(programVersion != null){
            out.append("\tVN:").append(programVersion);
        }
        if(commandLine != null){
            out.append("\tCL:").append(commandLine);
        }

    }
}
