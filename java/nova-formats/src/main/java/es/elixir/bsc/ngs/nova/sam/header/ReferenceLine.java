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

public class ReferenceLine extends AbstractHeaderLine {
    
    public final static String TAG = "@SQ";
    
    public final static String[] TAGS = 
                {"SN", "LN", "AH", "AN", "AS", "DS", "M5", "SP", "TP", "UR"};
    
    public final String referenceSequenceName;
    public final Long referenceLength;
    public final ReferenceLineAlternativeLocus alternativeLocus;
    public final String[] alternativeReferenceNames;
    public final String genomeAssemblyIdentifier;
    public final String description;
    public final String md5;
    public final String species;
    public final ReferenceLineMoleculeTopology moleculeTopology;
    public final String uri;

    public ReferenceLine(
            final String referenceSequenceName,
            final Long referenceLength,
            final ReferenceLineAlternativeLocus alternativeLocus,
            final String[] alternativeReferenceNames,
            final String genomeAssemblyIdentifier,
            final String description,
            final String md5,
            final String species,
            final ReferenceLineMoleculeTopology moleculeTopology,
            final String uri) {
        
        this.referenceSequenceName = referenceSequenceName;
        this.referenceLength = referenceLength;
        this.alternativeLocus = alternativeLocus;
        this.alternativeReferenceNames = alternativeReferenceNames;
        this.genomeAssemblyIdentifier = genomeAssemblyIdentifier;
        this.description = description;
        this.md5 = md5;
        this.species = species;
        this.moleculeTopology = moleculeTopology;
        this.uri = uri;
    }
    
    public ReferenceLine(final String line) {

        final String[] tags = parseHeaderLine(line, Arrays.copyOf(TAGS, TAGS.length));

        referenceSequenceName = tags[0];
        referenceLength = tags[1] == null ? null : Long.valueOf(tags[1]);
        alternativeLocus = tags[2] == null ? null : new ReferenceLineAlternativeLocus(tags[2]);
        alternativeReferenceNames = tags[3] == null ? null : tags[3].split(",");
        genomeAssemblyIdentifier = tags[4];
        description = tags[5];
        md5 = tags[6];
        species = tags[7];
        moleculeTopology = tags[8] == null ? null : ReferenceLineMoleculeTopology.valueOf(tags[8]);
        uri = tags[9];
    }

    @Override
    public void write(final PrintStream out) throws IOException {
        out.append(TAG);
        
        if (referenceSequenceName != null) {
            out.append("\tSN:").append(referenceSequenceName);
        }
        if(referenceLength != null) {
            out.append("\tLN:").print(referenceLength);
        }
        if(uri != null) {
            out.append("\tUR:").append(uri);
        }
        if(alternativeLocus != null) {
            out.append("\tAH:").print(alternativeLocus);
        }
        if(alternativeReferenceNames != null && alternativeReferenceNames.length > 0) {
            out.append("\tAN:").append(alternativeReferenceNames[0]);
            for(int i = 1; i < alternativeReferenceNames.length; i++) {
                out.append(",").append(alternativeReferenceNames[i]);
            }
        }
        if(genomeAssemblyIdentifier != null) {
            out.append("\tAS:").append(genomeAssemblyIdentifier);
        }
        if(description != null) {
            out.append("\tDS:").append(description);
        }
        if(md5 != null) {
            out.append("\tM5:").append(md5);
        }
        if(species != null) {
            out.append("\tSP:").append(species);
        }
        if(moleculeTopology != null) {
            out.append("\tTP:").print(moleculeTopology);
        }
    }
}
