package org.fhcrc.matsen.startree;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: cmccoy
 * Date: 12/3/13
 * Time: 5:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class SAMBEASTUtils {
    private SAMBEASTUtils() { throw new UnsupportedOperationException(); }

    private static class AlignedPair {
        private final int queryPosition;
        private final int referencePosition;

        private AlignedPair(int queryPosition, int referencePosition) {
            this.queryPosition = queryPosition;
            this.referencePosition = referencePosition;
        }

        private int getQueryPosition() {
            return queryPosition;
        }

        private int getReferencePosition() {
            return referencePosition;
        }

        public boolean consumesQuery() { return queryPosition >= 0; }
        public boolean consumesReference() { return referencePosition >= 0; }
    }

    private static List<AlignedPair> getAlignedPairs(final SAMRecord record) {
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

    /**
     * Create a BEAST-compatible alignment from a SAM record and reference bases.
     * @param record
     * @param rseq All reference sequence bases
     * @return
     */
    public static Alignment alignmentOfRecord(final SAMRecord record, final byte[] rseq) {
        StringBuilder raln = new StringBuilder(), qaln = new StringBuilder();
        SimpleAlignment alignment = new SimpleAlignment();

        final byte[] qseq = record.getReadBases();

        for(int i = 0; i < record.getAlignmentStart() - 1; i++) {
            raln.append((char)rseq[i]);
            qaln.append('-');
        }
        for(AlignedPair p : getAlignedPairs(record)) {
            if(p.consumesReference()) {
                raln.append((char)rseq[p.getReferencePosition()]);
                if(p.consumesQuery())
                    qaln.append((char)qseq[p.getQueryPosition()]);
                else
                    qaln.append('-');
            }
        }
        for(int i = qaln.length(); i < rseq.length; i++) {
            raln.append((char)rseq[i]);
            qaln.append('-');
        }

        alignment.addSequence(new Sequence(new Taxon(record.getReferenceName()), raln.toString()));
        alignment.addSequence(new Sequence(new Taxon(record.getReadName()), qaln.toString()));

        return alignment;
    }
}
