package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.alignmentsToSamTag;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_CIGAR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.markdups.common.PartitionData;
import com.hartwig.hmftools.markdups.umi.ConsensusReads;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MarkDuplicatesTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final RecordWriter mWriter;

    private final ChromosomeReader mChromosomeReader;
    private final ChromosomeReader mChrReaderUMIs;

    public MarkDuplicatesTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();

        MarkDupsConfig config = new MarkDupsConfig(1000, 1000, mRefGenome, false);
        mWriter = new RecordWriter(config);
        mWriter.setCacheReads();

        ChrBaseRegion chrBaseRegion = new ChrBaseRegion(CHR_1, 1, 100000);
        mChromosomeReader = new ChromosomeReader(chrBaseRegion, config, mWriter, new PartitionDataStore(config));

        MarkDupsConfig umiConfig = new MarkDupsConfig(1000, 1000, mRefGenome, true);
        mChrReaderUMIs = new ChromosomeReader(chrBaseRegion, umiConfig, mWriter, new PartitionDataStore(umiConfig));
    }

    @Test
    public void testNonDuplicateFragment()
    {
        // test non-duplicate fragments processing reads in various orders
        int readPos = 100;
        int matePos = 200;
        int suppPos = 2000;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setMateNegativeStrandFlag(true);
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mChromosomeReader.processRead(read1);
        mChromosomeReader.flushReadPositions();

        PartitionData partitionData = mChromosomeReader.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(1, partitionData.resolvedFragmentStateMap().size());
        assertEquals(1, mWriter.recordWriteCount());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        mChromosomeReader.processRead(mate1);
        assertEquals(1, partitionData.resolvedFragmentStateMap().size());
        assertEquals(2, mWriter.recordWriteCount());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChromosomeReader.processRead(supp1);
        mChromosomeReader.onChromosomeComplete();
        assertEquals(3, mWriter.recordWriteCount());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
    }

    @Test
    public void testSupplementaryFragments()
    {
        // first a supplementary linked to a distant mate but in the same partition as the primary
        int readPos = 100;
        int matePos = 2000;
        int suppPos = 200;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                false, null);

        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mChromosomeReader.processRead(read1);

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChromosomeReader.processRead(supp1);

        PartitionData partitionData = mChromosomeReader.partitionDataStore().getOrCreatePartitionData("1_0");
        partitionData.incompleteFragmentMap().containsKey(supp1.getReadName());

        mChromosomeReader.flushReadPositions();

        partitionData.incompleteFragmentMap().isEmpty();

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        mChromosomeReader.processRead(mate1);
        mChromosomeReader.onChromosomeComplete();

        assertEquals(3, mWriter.recordWriteCount());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
    }

        @Test
    public void testCandidateDuplicates()
    {
        // no use of mate CIGAR so rely on putting reads together to determine duplicate status
        int suppPos = 500;
        int readPos = 1500;
        int matePos = 2500;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setMateNegativeStrandFlag(true);

        SAMRecord read2 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setMateNegativeStrandFlag(true);

        // send through supplementaries first
        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChromosomeReader.processRead(supp1);
        mChromosomeReader.processRead(supp2);
        mChromosomeReader.flushPendingIncompletes();

        PartitionData partitionData = mChromosomeReader.partitionDataStore().getOrCreatePartitionData("1_1");

        assertEquals(2, partitionData.incompleteFragmentMap().size());
        assertEquals(0, mWriter.recordWriteCount());

        // now the primaries, as candidates
        mChromosomeReader.processRead(read1);
        mChromosomeReader.processRead(read2);
        mChromosomeReader.flushReadPositions();

        assertEquals(1, partitionData.candidateDuplicatesMap().size());
        assertEquals(2, partitionData.incompleteFragmentMap().size());
        assertEquals(0, mWriter.recordWriteCount());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate2.setFirstOfPairFlag(false);
        mate2.setSecondOfPairFlag(true);

        mChromosomeReader.processRead(mate1);
        assertEquals(0, mWriter.recordWriteCount());

        mChromosomeReader.processRead(mate2);
        mChromosomeReader.flushPendingIncompletes();
        assertEquals(6, mWriter.recordWriteCount());

        assertTrue(partitionData.candidateDuplicatesMap().isEmpty());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
    }

    @Test
    public void testUmiGroup()
    {
        // 2 primaries in a UMI group, followed by their mates in the same partition and then their supps in a different partition
        String umidId = "TATTAT";
        int readPos = 100;
        int matePos = 200;
        int suppPos = 2000;

        SAMRecord read1 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setMateNegativeStrandFlag(true);
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setMateNegativeStrandFlag(true);
        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.umiGroupMap().size());
        assertEquals(2, mWriter.recordWriteCount());
        assertEquals(1, mWriter.recordWriteCountConsensus());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate2.setFirstOfPairFlag(false);
        mate2.setSecondOfPairFlag(true);

        mChrReaderUMIs.processRead(mate1);
        assertEquals(2, mWriter.recordWriteCount());
        assertEquals(1, mWriter.recordWriteCountConsensus());

        mChrReaderUMIs.processRead(mate2);
        assertEquals(4, mWriter.recordWriteCount());
        assertEquals(2, mWriter.recordWriteCountConsensus());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(supp1);
        assertEquals(4, mWriter.recordWriteCount());
        assertEquals(2, mWriter.recordWriteCountConsensus());

        mChrReaderUMIs.processRead(supp2);
        mChrReaderUMIs.onChromosomeComplete();
        assertEquals(6, mWriter.recordWriteCount());
        assertEquals(3, mWriter.recordWriteCountConsensus());
        assertTrue(partitionData.umiGroupMap().isEmpty());
    }

    @Test
    public void testUmiGroupMultipleSupplementaries()
    {
        // one primary has multiple supplementaries and the other has a different supp
        String umidId = "TATTAT";
        int readPos = 100;
        int matePos = 200;
        int suppPos = 2000;

        SAMRecord read1 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null);

        List<SupplementaryReadData> read1Supps = Lists.newArrayList(
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                new SupplementaryReadData(CHR_2, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                new SupplementaryReadData(CHR_3, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setAttribute(SUPPLEMENTARY_ATTRIBUTE, alignmentsToSamTag(read1Supps));

        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos + 2, TEST_READ_BASES, "2S98M", CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_2, suppPos + 1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.umiGroupMap().size());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos + 2, true,
                false, null);
        mate2.setFirstOfPairFlag(false);
        mate2.setSecondOfPairFlag(true);

        mChrReaderUMIs.processRead(mate1);
        mChrReaderUMIs.processRead(mate2);

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read1Supp1);
        mChrReaderUMIs.onChromosomeComplete();

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_2, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord read2Supp1 = createSamRecord(
                read2.getReadName(), CHR_2, suppPos + 1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos + 2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read1Supp2);
        mChrReaderUMIs.processRead(read2Supp1);
        mChrReaderUMIs.onChromosomeComplete();

        SAMRecord read1Supp3 = createSamRecord(
                read1.getReadName(), CHR_3, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(read1Supp3);
        mChrReaderUMIs.onChromosomeComplete();

        assertEquals(7, mWriter.recordWriteCount());
        assertEquals(3, mWriter.recordWriteCountConsensus());

        partitionData.writeRemainingReads(mWriter, mChrReaderUMIs.consensusReads(), false);
        assertEquals(8, mWriter.recordWriteCount());
    }

    private String nextReadId(final String umiId)
    {
        return format("%s:%s", mReadIdGen.nextId(), umiId);

    }

}