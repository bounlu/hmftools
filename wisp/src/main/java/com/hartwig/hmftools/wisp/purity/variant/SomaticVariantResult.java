package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.util.StringJoiner;

import com.hartwig.hmftools.wisp.purity.ResultsWriter;

public class SomaticVariantResult
{
    public final int TotalVariants;
    public final int CalcVariants; // used in purity fit - passes filters and either avg qual > threshold or has no allele fragments
    public final int Frag1Variants;
    public final int Frag2PlusVariants;
    public final ClonalityMethod Method;
    public final int ClonalVariants;
    public final double ClonalDropoutRate;
    public final double WeightedAvgDepth;
    public final int TotalFragments;
    public final int UmiRefNonDual;
    public final int UmiRefDual;
    public final int AlleleFragments; //
    public final int UmiAlleleNonDual;
    public final int UmiAlleleDual;
    public final double DepthMedian;
    public final int NonZeroDepthCount;
    public final double NonZeroDepthMedian;

    public final double QualPerAdTotal;

    public final double TumorVaf;
    public final double AdjustedTumorVaf;
    public final double RawSomaticPurity;
    public final FragmentCalcResult AllFragsResult;
    public final FragmentCalcResult DualFragsResult;
    public final FragmentCalcResult LimitOfDetectionResult;

    private boolean mValid;

    public static final SomaticVariantResult INVALID_RESULT = new SomaticVariantResult(false);

    public SomaticVariantResult(
            boolean valid, int totalVariants, int calcVariants, int frag1Variants, int frag2PlusVariants, int clonalVariants,
            final ClonalityMethod method, final double dropoutRate, double weightedAvgDepth,
            final SomaticVariantCounts sampleCounts, final UmiTypeCounts umiTypeCounts,
            double qualPerAdTotal, double tumorVaf, double adjustedTumorVaf, double rawSomaticPurity,
            final FragmentCalcResult allFragsResult, final FragmentCalcResult dualFragsResult, final FragmentCalcResult lodFragsResult)
    {
        TotalVariants = totalVariants;
        CalcVariants = calcVariants;
        Frag1Variants = frag1Variants;
        Frag2PlusVariants = frag2PlusVariants;
        ClonalVariants = clonalVariants;
        ClonalDropoutRate = dropoutRate;
        Method = method;
        WeightedAvgDepth = weightedAvgDepth;
        TotalFragments = sampleCounts.totalFragments();
        UmiRefNonDual = umiTypeCounts.RefNone + umiTypeCounts.RefSingle;
        UmiRefDual = umiTypeCounts.RefDual;
        AlleleFragments = sampleCounts.alleleFragments();
        UmiAlleleNonDual = umiTypeCounts.AlleleNone + umiTypeCounts.AlleleSingle;
        UmiAlleleDual = umiTypeCounts.AlleleDual;
        QualPerAdTotal = qualPerAdTotal;
        DepthMedian = sampleCounts.medianDepth(false);
        NonZeroDepthCount = sampleCounts.nonZeroDepth();
        NonZeroDepthMedian = sampleCounts.medianDepth(true);
        TumorVaf = tumorVaf;
        AdjustedTumorVaf = adjustedTumorVaf;
        RawSomaticPurity = rawSomaticPurity;
        AllFragsResult = allFragsResult;
        DualFragsResult = dualFragsResult;
        LimitOfDetectionResult = lodFragsResult;
        mValid = valid;
    }

    public SomaticVariantResult(final boolean valid)
    {
        mValid = valid;
        TotalVariants = 0;
        CalcVariants = 0;
        Frag1Variants = 0;
        Frag2PlusVariants = 0;
        ClonalVariants = 0;
        Method = ClonalityMethod.NONE;
        ClonalDropoutRate = 0;
        WeightedAvgDepth = 0;
        TotalFragments = 0;
        UmiRefNonDual = 0;
        UmiRefDual = 0;
        UmiAlleleNonDual = 0;
        UmiAlleleDual = 0;
        AlleleFragments = 0;
        QualPerAdTotal = 0;
        DepthMedian = 0;
        NonZeroDepthCount = 0;
        NonZeroDepthMedian = 0;
        TumorVaf = 0;
        AdjustedTumorVaf = 0;
        RawSomaticPurity = 0;
        AllFragsResult = FragmentCalcResult.INVALID;
        DualFragsResult = FragmentCalcResult.INVALID;
        LimitOfDetectionResult = FragmentCalcResult.INVALID;
    }

    public boolean valid() { return mValid; }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("TotalVariants");
        sj.add("CalcVariants");
        sj.add("Frag1Variants");
        sj.add("Frag2PlusVariants");
        sj.add("ClonalPeakVariants");
        sj.add("ClonalMethod");
        sj.add("ClonalDropoutRate");
        sj.add("RawSomaticPurity");
        sj.add(FragmentCalcResult.header(""));
        sj.add(FragmentCalcResult.header("Dual"));
        sj.add("LodPurity");
        sj.add("TotalFragments");
        sj.add("UmiRefNonDual");
        sj.add("UmiRefDual");
        sj.add("AlleleFragments");
        sj.add("UmiAlleleNonDual");
        sj.add("UmiAlleleDual");
        sj.add("QualPerAdTotal");
        sj.add("TumorVaf");
        sj.add("AdjustedTumorVaf");
        sj.add("WeightedAvgDepth");
        sj.add("DepthMedian");
        sj.add("NonZeroDepthCount");
        sj.add("NonZeroDepthMedian");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(format("%d", TotalVariants));
        sj.add(format("%d", CalcVariants));
        sj.add(format("%d", Frag1Variants));
        sj.add(format("%d", Frag2PlusVariants));
        sj.add(format("%d", ClonalVariants));
        sj.add(String.valueOf(Method));
        sj.add(format("%.2f", ClonalDropoutRate));
        sj.add(ResultsWriter.formatPurityValue(RawSomaticPurity));
        sj.add(AllFragsResult.toTsv());
        sj.add(DualFragsResult.toTsv());
        sj.add(ResultsWriter.formatPurityValue(LimitOfDetectionResult.EstimatedPurity));

        // inputs
        sj.add(format("%d", TotalFragments));
        sj.add(format("%d", UmiRefNonDual));
        sj.add(format("%d", UmiRefDual));
        sj.add(format("%d", AlleleFragments));
        sj.add(format("%d", UmiAlleleNonDual));
        sj.add(format("%d", UmiAlleleDual));
        sj.add(format("%.1f", QualPerAdTotal));
        sj.add(ResultsWriter.formatPurityValue(TumorVaf));
        sj.add(ResultsWriter.formatPurityValue(AdjustedTumorVaf));
        sj.add(format("%.0f", WeightedAvgDepth));
        sj.add(format("%.1f", DepthMedian));
        sj.add(format("%d", NonZeroDepthCount));
        sj.add(format("%.1f", NonZeroDepthMedian));

        return sj.toString();
    }
}