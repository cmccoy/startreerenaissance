package org.fhcrc.matsen.startree;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by cmccoy on 12/5/13.
 */
public class SAMUtils {
    public static List<AlignedPair> getAlignedPairs(final SAMRecord record) {
        int q = 0, r = record.getAlignmentStart() - 1;
        List<AlignedPair> result = new ArrayList<AlignedPair>();
        for(final CigarElement e : record.getCigar().getCigarElements()) {
            final CigarOperator op = e.getOperator();
            for(int i = 0; i < e.getLength(); i++) {
                result.add(new AlignedPair(op.consumesReadBases() ? q++ : -1,
                        op.consumesReferenceBases() ? r++ : -1));
            }
        }
        return result;
    }
}
