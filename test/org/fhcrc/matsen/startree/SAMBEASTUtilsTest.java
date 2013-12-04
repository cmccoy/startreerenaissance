package org.fhcrc.matsen.startree;

import dr.evolution.alignment.Alignment;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Test;

/**
 * Created with IntelliJ IDEA.
 * User: cmccoy
 * Date: 12/3/13
 * Time: 6:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class SAMBEASTUtilsTest {
    private final byte[] reference = "ACCGTACTA".getBytes();
    private final String referenceName = "Ref1";
    private final String queryName = "QSEQ1";
    private SAMFileHeader header;

    @Before
    public void setUp() throws Exception {
        header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord(referenceName, reference.length));
    }

    @Test
    public void testAlignmentOfRecord() throws Exception {
        SAMRecord record = new SAMRecord(this.header);
        record.setReadBases("TCCGTTAGTC".getBytes());
        record.setReadName(queryName);
        record.setCigarString("1S3M1I4M1S");
        record.setReferenceIndex(0);
        record.setReferenceName(this.referenceName);
        record.setAlignmentStart(2);

        assertEquals(referenceName, record.getReferenceName());

        Alignment result = SAMBEASTUtils.alignmentOfRecord(record, this.reference);
        assertEquals(2, result.getTaxonCount());
        assertEquals(new String(reference), result.getAlignedSequenceString(0));
        assertEquals("-CCGTAGT-", result.getAlignedSequenceString(1));
        assertEquals(referenceName, result.getTaxonId(0));
        assertEquals(queryName, result.getTaxonId(1));
    }
}

