package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.MAX_MAP_QUALITY;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;

public class QualityCalculator
{
    private final QualityConfig mConfig;
    private final BqrRecordMap mQualityRecalibrationMap;
    private final IndexedBases mRefBases;

    private static final int MAX_HIGHLY_POLYMORPHIC_GENES_QUALITY = 10;

    public QualityCalculator(
            final QualityConfig config, final BqrRecordMap qualityRecalibrationMap, final IndexedBases refBases)
    {
        mConfig = config;
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mRefBases = refBases;
    }

    public static int modifiedMapQuality(
            final QualityConfig config, final GenomePosition position, int mapQuality, double readEvents, boolean isImproperPair)
    {
        if(config.isHighlyPolymorphic(position))
        {
            return min(MAX_HIGHLY_POLYMORPHIC_GENES_QUALITY, mapQuality - config.FixedPenalty);
        }

        int improperPairPenalty = isImproperPair ? config.ImproperPairPenalty : 0;
        int eventPenalty = (int)round(max(0, readEvents - 1) * config.ReadEventsPenalty);

        int modifiedMapQuality = mapQuality - config.FixedPenalty - improperPairPenalty - eventPenalty;

        return config.MapQualityRatioFactor > 0 ? min(MAX_MAP_QUALITY, modifiedMapQuality) : modifiedMapQuality;
    }

    public static class QualityScores
    {
        public final int ModifiedMapQuality;
        public final double ModifiedBaseQuality;
        public final double ModifiedQuality;

        private double mRawBaseQuality;

        public QualityScores(int modifiedMapQuality, double modifiedBaseQuality, double modifiedQuality)
        {
            ModifiedMapQuality = modifiedMapQuality;
            ModifiedBaseQuality = modifiedBaseQuality;
            ModifiedQuality = modifiedQuality;

            mRawBaseQuality = 0.0;
        }

        public double rawBaseQuality()
        {
            return mRawBaseQuality;
        }

        public void setRawBaseQuality(double rawBaseQuality)
        {
            mRawBaseQuality = rawBaseQuality;
        }
    }

    public QualityScores calculateQualityScores(
            final ReadContextCounter readContextCounter, int readBaseIndex, final SAMRecord record, double numberOfEvents, double rawBaseQuality)
    {
        double baseQuality = readContextCounter.isIndel() ?
                rawBaseQuality : baseQuality(readContextCounter, readBaseIndex, record, readContextCounter.variant().ref().length());

        int mapQuality = record.getMappingQuality();
        boolean isImproperPair = isImproperPair(record);

        int modifiedMapQuality = modifiedMapQuality(mConfig, readContextCounter.variant(), mapQuality, numberOfEvents, isImproperPair);

        double modifiedBaseQuality = baseQuality - mConfig.BaseQualityFixedPenalty;

        if(mConfig.DistanceFromReadEdgeFactor > 0)
        {
            int distanceFromReadEdge = readDistanceFromEdge(readContextCounter, readBaseIndex, record);
            int readEdgePenalty = max(mConfig.DistanceFromReadEdgeFactor * distanceFromReadEdge - mConfig.DistanceFromReadEdgeFixedPenalty, 0);
            modifiedBaseQuality = min(modifiedBaseQuality, readEdgePenalty);
        }

        double modifiedQuality = max(0, min(modifiedMapQuality, modifiedBaseQuality));

        /*
        if(readContextCounter.logEvidence() && !SG_LOGGER.isTraceEnabled())
        {
            SG_LOGGER.trace(format("variant(%s) read(%s) distFromEdge(%d) events(%.1f) qual(map=%d rawBase=%.1f base=%.1f) modified(map=%d base=%.1f)",
                    readContextCounter.varString(), record.getReadName(), distanceFromReadEdge, numberOfEvents,
                    mapQuality, rawBaseQual, baseQuality, modifiedMapQuality, modifiedBaseQuality));
        }
        */

        return new QualityScores(max(0, modifiedMapQuality), max(0.0, modifiedBaseQuality), modifiedQuality);
    }

    public static boolean isImproperPair(final SAMRecord record) { return record.getReadPairedFlag() && !record.getProperPairFlag(); }

    public static double rawBaseQuality(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        if(readContextCounter.isIndel())
            return readContextCounter.readContext().avgCentreQuality(readIndex, record);

        if(readContextCounter.isSnv())
            return record.getBaseQualities()[readIndex];

        int varLength = readContextCounter.variant().ref().length();

        double baseQualTotal = 0;

        for(int i = readIndex; i < readIndex + varLength; ++i)
        {
            baseQualTotal += record.getBaseQualities()[i];
        }

        return baseQualTotal / varLength;
    }

    private double baseQuality(final ReadContextCounter readContextCounter, int startReadIndex, final SAMRecord record, int length)
    {
        if(readContextCounter.isSnv())
        {
            // simplified version of the MNV case below
            byte rawQuality = record.getBaseQualities()[startReadIndex];
            return recalibrateQuality(readContextCounter, readContextCounter.position(), 0, rawQuality);
        }

        // MNV case
        int maxIndex = min(startReadIndex + length, record.getBaseQualities().length) - 1;
        int maxLength = maxIndex - startReadIndex + 1;

        double quality = Integer.MAX_VALUE;
        for(int i = 0; i < maxLength; i++)
        {
            int refPosition = readContextCounter.position() + i;
            int readIndex = startReadIndex + i;
            byte rawQuality = record.getBaseQualities()[readIndex];

            double recalibratedQual = recalibrateQuality(readContextCounter, refPosition, i, rawQuality);
            quality = min(quality, recalibratedQual);
        }

        return quality;
    }

    private double recalibrateQuality(final ReadContextCounter readContextCounter, int refPosition, int refAltPos, byte rawQuality)
    {
        if(rawQuality == 0)
            return 0; // never adjust a zero qual up

        if(mQualityRecalibrationMap == null)
            return rawQuality;

        byte[] trinucleotideContext = mRefBases.trinucleotideContext(refPosition);

        return mQualityRecalibrationMap.getQualityAdjustment(
                (byte) readContextCounter.ref().charAt(refAltPos),
                (byte) readContextCounter.alt().charAt(refAltPos),
                trinucleotideContext, rawQuality);
    }

    private int readDistanceFromEdge(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        // calculate the left and right core positions in the context of this read
        int index = readContextCounter.readContext().readBasesPositionIndex();
        int leftIndex = readContextCounter.readContext().readBasesLeftCentreIndex();
        int rightIndex = readContextCounter.readContext().readBasesRightCentreIndex();

        int leftOffset = index - leftIndex;
        int rightOffset = rightIndex - index;

        int adjustedLeftIndex = readIndex - leftOffset;
        int adjustedRightIndex = readIndex + rightOffset;

        // take the smaller of the left and right core index
        return max(0, min(adjustedLeftIndex, record.getReadBases().length - 1 - adjustedRightIndex));
    }
}
