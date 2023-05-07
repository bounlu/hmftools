package com.hartwig.hmftools.common.test;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public final class SamRecordTestUtils
{
    public static final int DEFAULT_BASE_QUAL = 37;
    public static final String NO_MATE_CHR = "*";

    public static SAMSequenceDictionary SAM_DICTIONARY_V37;

    static
    {
        SAM_DICTIONARY_V37 = new SAMSequenceDictionary();

        RefGenomeCoordinates v37Coords = RefGenomeCoordinates.COORDS_37;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            SAM_DICTIONARY_V37.addSequence(new SAMSequenceRecord(chromosome.toString(), v37Coords.Lengths.get(chromosome)));
        }
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar, final String mateChr,
            int mateStart, boolean isReversed, boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.getHeader().setSequenceDictionary(SAM_DICTIONARY_V37);
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        HumanChromosome chromosome = HumanChromosome.fromString(chrStr);

        SAMRecord record = recordBuilder.addFrag(
                readId, chromosome.ordinal(), readStart, isReversed, false,
                cigar, readBases, DEFAULT_BASE_QUAL, false);

        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = DEFAULT_BASE_QUAL;

        record.setBaseQualities(qualities);
        record.setReferenceName(chrStr);
        record.setReferenceIndex(chromosome.ordinal()); // need to override since no header is present

        if(!mateChr.isEmpty() && !mateChr.equals(NO_MATE_CHR))
        {
            record.setMateReferenceName(mateChr);
            record.setMateAlignmentStart(mateStart);
            record.setMateReferenceIndex(HumanChromosome.fromString(mateChr).ordinal());
            record.setReadPairedFlag(true);
            record.setProperPairFlag(true);
        }
        else
        {
            record.setMateUnmappedFlag(true);
            record.setProperPairFlag(false);
        }

        // to be correct this should match the cigar element count
        record.setFirstOfPairFlag(true);

        record.setSupplementaryAlignmentFlag(isSupplementary);

        if(suppAlignment != null)
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, suppAlignment.asSamTag());

        if(chrStr.equals(mateChr))
            record.setInferredInsertSize(abs(readStart - mateStart));
        else
            record.setInferredInsertSize(0);

        return record;
    }

}