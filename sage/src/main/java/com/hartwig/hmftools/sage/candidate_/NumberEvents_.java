package com.hartwig.hmftools.sage.candidate_;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;

import com.hartwig.hmftools.sage.common.RefSequence;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Utility methods for calculating number of events in a read.
 */
public class NumberEvents_
{
    /**
     * Returns NM attribute adjusted so that indels count as one regardless of their length.
     */
    public static int calc(final SAMRecord record, final RefSequence refGenome)
    {
        int nm = rawNM(record, refGenome);

        int additionalIndels = 0;

        for(CigarElement cigarElement : record.getCigar())
        {
            switch(cigarElement.getOperator())
            {
                case D:
                case I:
                    additionalIndels += cigarElement.getLength() - 1;
                    break;
            }
        }

        return nm - additionalIndels;
    }

    /**
     * NM attribte of a read.
     */
    public static int rawNM(final SAMRecord record, final RefSequence refGenome)
    {
        Object nm = record.getAttribute(NUM_MUTATONS_ATTRIBUTE);
        if(nm instanceof Integer)
        {
            return (int) nm;
        }

        int offset = refGenome.alignment().Position - refGenome.alignment().Index - 1;
        return SequenceUtil.calculateSamNmTag(record, refGenome.alignment().Bases, offset);
    }

    /**
     * Return scaled a clipped number of soft clipped bases.
     */
    public static double calcSoftClipAdjustment(final SAMRecord record)
    {
        int softClippedBases = leftSoftClipLength(record) + rightSoftClipLength(record);

        return softClippedBases > 0 ? max(1, softClippedBases / SC_READ_EVENTS_FACTOR) : 0;
    }

    /**
     * Adjusts number of events for an MNV, so that the MVN only counts as a single event instead of one for each mutation in the MNV.
     */
    public static int calcWithMnvRaw(int numberOfEvents, final String ref, final String alt)
    {
        // Number of events includes each SNV as an additional event. This unfairly penalises MNVs.
        int differentBases = 0;
        for(int i = 0; i < alt.length(); i++)
        {
            if(alt.charAt(i) != ref.charAt(i))
            {
                differentBases++;
            }
        }

        // We subtract one later when we actually use this value so we need to add one back in here to be consistent with SNVs and INDELs
        return numberOfEvents - differentBases + 1;
    }
}
