package org.fhcrc.matsen.startree;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by cmccoy on 12/5/13.
 */
public class SAMUtils {
    /**
     * Get aligned pairs from a SAM record
     * <p/>
     * This mimics pysams AlignedRead.aligned_pairs
     *
     * @param record Aligned read
     * @return List of aligned (reference, query) indices
     */
    public static List<AlignedPair> getAlignedPairs(final SAMRecord record) {
        Preconditions.checkNotNull(record);
        int q = 0, r = record.getAlignmentStart() - 1;
        List<AlignedPair> result = Lists.newArrayList();
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
     * Read all sequences from a FASTA file, returning a map from name to sequence
     *
     * @param path Path to FASTA file
     * @return Map from name to sequence
     */
    public static Map<String, byte[]> readAllFasta(final File path) {
        Preconditions.checkNotNull(path);
        final Map<String, byte[]> result = new HashMap<String, byte[]>();

        final FastaSequenceFile file = new FastaSequenceFile(path, true);
        ReferenceSequence sequence = file.nextSequence();
        while (sequence != null) {
            result.put(sequence.getName(), sequence.getBases());
            sequence = file.nextSequence();
        }

        return result;
    }
}
