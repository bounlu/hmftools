package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class BaseCoverage
{
    private final MetricsConfig mConfig;
    private final int mRegionSize;
    private final int mRegionStart;
    private final List<ChrBaseRegion> mUnmappableRegions;

    private final int[] mBaseDepth;
    private final long[] mFilterTypeCounts;

    public BaseCoverage(final MetricsConfig config, int regionStart, int regionEnd, final List<ChrBaseRegion> unmappableRegions)
    {
        mConfig = config;
        mRegionSize = regionEnd - regionStart + 1;
        mRegionStart = regionStart;
        mBaseDepth = new int[mRegionSize];
        mFilterTypeCounts = new long[FilterType.values().length];
        mUnmappableRegions = unmappableRegions;
    }

    public void processRead(final SAMRecord read, final List<int[]> mateBaseCoords)
    {
        // some filters exclude all matched bases
        int alignedBases = read.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        boolean isConsensus = read.hasAttribute(CONSENSUS_READ_ATTRIBUTE);
        // the order in which the filters are applied matters and matches Picard CollectWgsMetrics

        if(read.getMappingQuality() < mConfig.MapQualityThreshold)
        {
            if(!isConsensus)
            {
                mFilterTypeCounts[FilterType.LOW_MAP_QUAL.ordinal()] += alignedBases;
            }
            return;
        }

        if(isConsensus)
        {
            mFilterTypeCounts[FilterType.DUPLICATE.ordinal()] -= alignedBases;
        }
        else if(read.getDuplicateReadFlag())
        {
            mFilterTypeCounts[FilterType.DUPLICATE.ordinal()] += alignedBases;
            return;
        }

        if(!read.getReadPairedFlag() || read.getMateUnmappedFlag())
        {
            mFilterTypeCounts[FilterType.MATE_UNMAPPED.ordinal()] += alignedBases;
            return;
        }

        int position = read.getAlignmentStart();
        int readIndex = 0;

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            switch(element.getOperator())
            {
                case S:
                    readIndex += element.getLength();
                    break;

                case D:
                case N:
                    position += element.getLength();
                    break;

                case H:
                    break;

                case I:
                    readIndex += element.getLength();
                    break;

                case M:
                    processMatchedBases(read, position, readIndex, element.getLength(), mateBaseCoords);

                    position += element.getLength();
                    readIndex += element.getLength();
                    break;

                default:
                    break;
            }
        }
    }

    private void processMatchedBases(
            final SAMRecord read, int posStart, int readIndexStart, int matchLength, final List<int[]> mateBaseCoords)
    {
        boolean checkUnmappable = !mUnmappableRegions.isEmpty()
                && mUnmappableRegions.stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), posStart, read.getAlignmentEnd()));

        for(int i = 0; i < matchLength; ++i)
        {
            int position = posStart + i;

            if(position < mRegionStart)
                continue;

            if(position >= mRegionStart + mRegionSize)
                break;

            if(checkUnmappable && mUnmappableRegions.stream().anyMatch(x -> x.containsPosition(position)))
                continue;

            int readIndex = readIndexStart + i;
            int baseIndex = position - mRegionStart;

            boolean lowBaseQual = read.getBaseQualities()[readIndex] < mConfig.BaseQualityThreshold;

            boolean overlapped = mateBaseCoords != null && mateBaseCoords.stream().anyMatch(x -> positionWithin(position, x[SE_START], x[SE_END]));

            boolean exceedsCoverage = mBaseDepth[baseIndex] >= mConfig.MaxCoverage;

            boolean passFilters = !lowBaseQual && !exceedsCoverage && !overlapped;

            if(passFilters)
            {
                ++mBaseDepth[baseIndex];
                ++mFilterTypeCounts[FilterType.UNFILTERED.ordinal()];
            }
            else
            {
                if(lowBaseQual)
                    ++mFilterTypeCounts[FilterType.LOW_BASE_QUAL.ordinal()];
                else if(overlapped)
                    ++mFilterTypeCounts[FilterType.OVERLAPPED.ordinal()];
                else if(exceedsCoverage)
                    ++mFilterTypeCounts[FilterType.MAX_COVERAGE.ordinal()];
            }
        }
    }

    public CoverageMetrics createMetrics()
    {
        CoverageMetrics metrics = new CoverageMetrics(mConfig.MaxCoverage);

        long coverageBases = 0;

        for(int i = 0; i < mBaseDepth.length; ++i)
        {
            int coverage = mBaseDepth[i];

            if(!mUnmappableRegions.isEmpty())
            {
                int position = mRegionStart + i;

                if(mUnmappableRegions.stream().anyMatch(x -> x.containsPosition(position)))
                    continue;
            }

            if(coverage == 0)
            {
                if(mConfig.ExcludeZeroCoverage)
                    continue;
            }
            else
            {
                ++coverageBases;
            }

            ++metrics.CoverageFrequency[coverage];
        }

        for(FilterType type : FilterType.values())
        {
            metrics.FilterTypeCounts[type.ordinal()] += mFilterTypeCounts[type.ordinal()];
        }

        metrics.addCoverageBases(coverageBases);

        return metrics;
    }

    public void clear()
    {
        for(int i = 0; i < mBaseDepth.length; ++i)
        {
            mBaseDepth[i] = 0;
        }

        for(int i = 0; i < mFilterTypeCounts.length; ++i)
        {
            mFilterTypeCounts[i] = 0;
        }
    }

    @VisibleForTesting
    public int[] baseDepth() { return mBaseDepth; }
}
