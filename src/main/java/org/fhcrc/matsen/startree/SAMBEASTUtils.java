package org.fhcrc.matsen.startree;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
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
    private SAMBEASTUtils() {
        throw new UnsupportedOperationException();
    }

    private static List<AlignedPair> getAlignedPairs(final SAMRecord record) {
        int q = 0, r = record.getAlignmentStart() - 1;
        List<AlignedPair> result = new ArrayList<AlignedPair>();
        for (final CigarElement e : record.getCigar().getCigarElements()) {
            final CigarOperator op = e.getOperator();
            for (int i = 0; i < e.getLength(); i++) {
                result.add(new AlignedPair(op.consumesReadBases() ? q++ : -1,
                        op.consumesReferenceBases() ? r++ : -1));
            }
        }
        return result;
    }


    /**
     * Create a BEAST-compatible alignment from a SAM record and reference bases.
     *
     * @param record Aligned read
     * @param rSeq   All reference sequence bases
     * @param onlyCertain Convert uncertain bases to 'N'?
     * @return An alignment, covering all bases of rSeq
     */
    public static Alignment alignmentOfRecord(final SAMRecord record, final byte[] rSeq, final boolean onlyCertain) {
        StringBuilder rAln = new StringBuilder(), qAln = new StringBuilder();
        SimpleAlignment alignment = new SimpleAlignment();

        final byte[] qSeq = record.getReadBases();
        byte[] certainty = null;
        if (onlyCertain)
            certainty = record.getByteArrayAttribute("bq");


        for (int i = 0; i < record.getAlignmentStart() - 1; i++) {
            rAln.append((char) rSeq[i]);
            qAln.append('-');
        }
        for (AlignedPair p : getAlignedPairs(record)) {
            if (p.consumesReference()) {
                rAln.append((char) rSeq[p.getReferencePosition()]);
                if (p.consumesQuery()) {
                    if(onlyCertain && certainty[p.getQueryPosition()] % 100 != 0)
                        qAln.append('N');
                    else
                        qAln.append((char) qSeq[p.getQueryPosition()]);
                } else
                    qAln.append('-');
            }
        }
        for (int i = qAln.length(); i < rSeq.length; i++) {
            rAln.append((char) rSeq[i]);
            qAln.append('-');
        }

        alignment.addSequence(new Sequence(new Taxon(record.getReferenceName()), rAln.toString()));
        alignment.addSequence(new Sequence(new Taxon(record.getReadName()), qAln.toString()));

        return alignment;
    }

    /**
     * Create a BEAST-compatible alignment from a SAM record and reference bases.
     *
     * @param record Aligned read
     * @param rSeq   All reference sequence bases
     * @return An alignment, covering all bases of rSeq
     */
    public static Alignment alignmentOfRecord(final SAMRecord record, final byte[] rSeq) {
        return alignmentOfRecord(record, rSeq, false);
    }
}
