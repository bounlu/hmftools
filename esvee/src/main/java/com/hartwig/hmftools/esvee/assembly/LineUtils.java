package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_TEST_LEN;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.getReadIndexAtReferencePosition;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.Junction;

public final class LineUtils
{
    public static boolean isLineSequence(final byte[] bases, final int indexStart, final int indexEnd)
    {
        // returns true if the whole sequence matches a LINE base
        if(indexStart < 0 || indexEnd >= bases.length)
            return false;

        byte lineBase = bases[indexStart];

        if(lineBase != LINE_BASE_A && lineBase != LINE_BASE_T)
            return false;

        for(int i = indexStart + 1; i <= indexEnd; ++i)
        {
            if(bases[i] != lineBase)
                return false;
        }

        return true;
    }

    public static int findLineSequenceCount(final byte[] bases, final int indexStart, final int indexEnd, final byte lineBase)
    {
        if(indexStart < 0 || indexEnd >= bases.length)
            return 0;

        if(indexEnd - indexStart + 1 < LINE_POLY_AT_REQ)
            return 0;

        int lineBaseCount = 0;
        int otherCount = 0;
        for(int i = indexStart; i <= indexEnd; ++i)
        {
            if(bases[i] == lineBase)
            {
                ++lineBaseCount;
            }
            else
            {
                ++otherCount;

                if(otherCount > LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ)
                    break;
            }
        }

        return lineBaseCount >= LINE_POLY_AT_REQ ? lineBaseCount : 0;
    }

    public static int findConsensusLineExtension(final List<Read> reads, final Junction junction)
    {
        // find the poly A/T sequence by first favouring reads heading to the 5' end, ending in non-poly-T bases, taking their median bounds
        Map<Integer,Integer> lengthFrequency = findLineExtensionFrequency(reads, junction, true, true);

        if(!lengthFrequency.isEmpty())
            return calcMedian(lengthFrequency);

        lengthFrequency = findLineExtensionFrequency(reads, junction, false, true);

        if(!lengthFrequency.isEmpty())
            return calcMedian(lengthFrequency);

        lengthFrequency = findLineExtensionFrequency(reads, junction, false, false);

        return calcMedian(lengthFrequency);
    }

    private static final int MAX_NON_LINE_BASES = LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ;

    @VisibleForTesting
    public static Map<Integer,Integer> findLineExtensionFrequency(
            final List<Read> reads, final Junction junction, boolean fivePrimeOnly, boolean terminatingOnly)
    {
        byte lineBase = junction.isForward() ? LINE_BASE_T : LINE_BASE_A;
        boolean isForward = junction.isForward();
        int maxNonLine = LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ;

        Map<Integer,Integer> lengthFrequency = Maps.newHashMap();

        for(Read read : reads)
        {
            if(!read.hasLineTail())
                continue;

            if(fivePrimeOnly)
            {
                // only count if moving towards the 5' end of the read
                if(read.positiveStrand() == isForward)
                    continue;
            }

            int readJunctionIndex = getReadIndexAtReferencePosition(read, junction.Position, true);

            if(readJunctionIndex == INVALID_INDEX)
                continue;

            int firstExtensionIndex = readJunctionIndex + (isForward ? 1 : -1);
            int lineExtensionIndex = findLineExtensionEndIndex(read, lineBase, firstExtensionIndex, isForward);
            int lineLength = abs(firstExtensionIndex - lineExtensionIndex) + 1;

            if(lineLength < LINE_POLY_AT_REQ)
                continue;

            if(terminatingOnly)
            {
                // compare LINE length vs soft-clip length
                int softClipLength = isForward ? read.rightClipLength() : read.leftClipLength();

                if(softClipLength == lineLength)
                    continue;
            }

            Integer count = lengthFrequency.get(lineLength);
            lengthFrequency.put(lineLength, count != null ? count + 1 : 1);
        }

        return lengthFrequency;
    }

    public static int findLineExtensionEndIndex(final Read read, final byte lineBase, int extensionIndex, boolean isForward)
    {
        int index = extensionIndex;
        int lineLength = 0;
        int nonLineCount = 0;

        while(index >= 0 && index < read.basesLength())
        {
            if(read.getBases()[index] == lineBase)
            {
                ++lineLength;
            }
            else
            {
                // break on any mismatch once the minimum required LINE site length has been reached
                if(lineLength >= LINE_POLY_AT_REQ)
                    break;

                if(aboveMinQual(read.getBaseQuality()[index]))
                {
                    ++nonLineCount;

                    if(nonLineCount > MAX_NON_LINE_BASES)
                        break;
                }
            }

            index += isForward ? 1 : -1;
        }

        // back up to end of read or last line base
        index += isForward ? -1 : 1;
        return index;
    }

    private static int calcMedian(final Map<Integer,Integer> lengthFrequency)
    {
        if(lengthFrequency.isEmpty())
            return 0;

        int totalEntries = lengthFrequency.values().stream().mapToInt(x -> x.intValue()).sum();

        List<Integer> lengths = Lists.newArrayListWithCapacity(totalEntries);

        for(Map.Entry<Integer,Integer> entry : lengthFrequency.entrySet())
        {
            for(int i = 0; i < entry.getValue(); ++i)
            {
                lengths.add(entry.getKey());
            }
        }

        Collections.sort(lengths);

        int medianIndex = totalEntries / 2;
        return lengths.get(medianIndex);
    }
}