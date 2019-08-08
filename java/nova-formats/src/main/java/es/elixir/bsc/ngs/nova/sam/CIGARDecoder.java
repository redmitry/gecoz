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

import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Dmitry Repchevsky
 */

public class CIGARDecoder {
    private final static Pattern PATTERN = Pattern.compile("([0-9]+)([MIDNSHP=X]{1})");

    private final byte[] operation;
    private final int[] operationLength;

    public CIGARDecoder(final String cigar) {
        byte[] operationTmp = new byte[1];
        int[] operationLengthTmp = new int[1];
        final Matcher m = PATTERN.matcher(cigar);

        while (m.find()) {
            final String ln = m.group(1);
            final String ch = m.group(2);
            if (ln != null && ch != null) {
                operationTmp = Arrays.copyOf(operationTmp, operationTmp.length + 1);
                operationLengthTmp = Arrays.copyOf(operationLengthTmp, operationLengthTmp.length + 1);


                operationTmp[operationTmp.length - 1] = ch.getBytes()[0];
                operationLengthTmp[operationTmp.length - 1] = Integer.parseInt(ln);
            }
        }
        operation = operationTmp;
        operationLength = operationLengthTmp;
    }

    public final byte[] getOperation() {
        return operation;
    }

    public final int[] getOperationLength() {
        return operationLength;
    }
}