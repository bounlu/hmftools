package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.utils.pcf.PCFSource.TUMOR_BAF;
import static com.hartwig.hmftools.common.utils.pcf.PCFSource.TUMOR_RATIO;
import static com.hartwig.hmftools.purple.PurpleTestUtils.REF_SAMPLE_ID;
import static com.hartwig.hmftools.purple.PurpleTestUtils.TUMOR_SAMPLE_ID;
import static com.hartwig.hmftools.purple.PurpleTestUtils.buildCobaltChromosomes;
import static com.hartwig.hmftools.purple.PurpleTestUtils.buildDefaultConfigBuilder;
import static com.hartwig.hmftools.purple.PurpleTestUtils.createAmberBaf;
import static com.hartwig.hmftools.purple.PurpleTestUtils.createCobaltRatio;
import static com.hartwig.hmftools.purple.PurpleTestUtils.createObservedRegion;
import static com.hartwig.hmftools.purple.PurpleTestUtils.createSegmentation;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.FittingConfig.MAX_PLOIDY;
import static com.hartwig.hmftools.purple.config.FittingConfig.PURITY_INCREMENT;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.purple.PurpleTestUtils;
import com.hartwig.hmftools.purple.config.AmberData;
import com.hartwig.hmftools.purple.config.CobaltData;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.config.SampleData;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.Segmentation;
import com.hartwig.hmftools.purple.somatic.SomaticVariantCache;
import com.hartwig.hmftools.purple.sv.SomaticSvCache;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurityPloidyFitTest
{
    private final PurpleConfig mConfig;
    private final SampleData mSampleData;
    private final AmberData mAmberData;
    private final CobaltData mCobaltData;

    private final SomaticSvCache mSvCache;
    private final SomaticVariantCache mSomaticCache;

    private final RegionFitCalculator mRegionFitCalculator;

    private final ReferenceData mReferenceData;
    private final Segmentation mSegmentation;

    public PurityPloidyFitTest()
    {
        ConfigBuilder configBuilder = buildDefaultConfigBuilder();

        configBuilder.setValue(PURITY_INCREMENT, 0.2);
        configBuilder.setValue(MAX_PLOIDY, 4);
        configBuilder.setValue(PURITY_INCREMENT, 0.2);

        mConfig = PurpleTestUtils.buildPurpleConfig(configBuilder);

        mAmberData = new AmberData(100, Gender.MALE);

        CobaltChromosomes cobaltChromosomes = buildCobaltChromosomes();
        mCobaltData = new CobaltData(cobaltChromosomes);

        mSomaticCache = new SomaticVariantCache(mConfig);
        mSvCache = new SomaticSvCache();

        mRegionFitCalculator = new RegionFitCalculator(cobaltChromosomes, mConfig.Fitting, mAmberData.AverageTumorDepth);

        mReferenceData = new ReferenceData(mConfig);
        mSegmentation = createSegmentation(mReferenceData);

        mSampleData = new SampleData(REF_SAMPLE_ID, TUMOR_SAMPLE_ID, mAmberData, mCobaltData, mSvCache, mSomaticCache);

        buildDefaultProfile();
    }

    private void buildDefaultProfile()
    {
        Chromosome chr1 = HumanChromosome._1;

        // Amber data
        mAmberData.ChromosomeBafs.put(chr1, createAmberBaf(chr1.toString(), 1000, 1.0, 0.5));

        mAmberData.TumorSegments.put(chr1, new PCFPosition(TUMOR_BAF, chr1.toString(), 100000));

        // Cobalt data
        mCobaltData.TumorSegments.put(chr1, new PCFPosition(TUMOR_RATIO, chr1.toString(), 100000));

        List<CobaltRatio> cobaltRatios = Lists.newArrayList();
        cobaltRatios.add(createCobaltRatio(chr1.toString(), 1000, 0.5));
        cobaltRatios.add(createCobaltRatio(chr1.toString(), 2000, 0.5));
        cobaltRatios.add(createCobaltRatio(chr1.toString(), 3000, 1));
        cobaltRatios.add(createCobaltRatio(chr1.toString(), 4000, 0.5));
        mCobaltData.Ratios.put(chr1, cobaltRatios);
    }

    
    // inspired by NO_TUMOR sample WIDE01010763T
    private List<ObservedRegion> buildNoTumorObservedRegions()
    {
        List<ObservedRegion> observedRegions = Lists.newArrayList();

        Chromosome chr1 = HumanChromosome._1;
        
        observedRegions.add(createObservedRegion(
                chr1.toString(), 1, 1000, 0.4466, 0.4138));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 1001, 2000, 0.4655, 0.5417));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 2001, 3000, 0.4316, 0.4783));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 3001, 4000, 0.4211, 0.4348));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 4001, 5000, 0.3830, 0.4667));

        return observedRegions;
    }

    // inspired by NORMAL (copyNumber mode) sample COREDB011180T
    private List<ObservedRegion> buildCopyNumberObservedRegions()
    {
        List<ObservedRegion> observedRegions = Lists.newArrayList();

        Chromosome chr1 = HumanChromosome._1;

        observedRegions.add(createObservedRegion(
                chr1.toString(), 1, 1000, 0.4286, 0.4643));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 1001, 2000, 0.4699, 0.5952));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 2001, 3000, 0.4787, 0.4524));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 3001, 4000, 0.4725, 0.6111));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 4001, 5000, 0.4937, 0.6250));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 5001, 6000, 0.4348, 0.5588));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 6001, 7000, 0.5275, 0.5152));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 7001, 8000, 0.500, 0.5294));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 8001, 9000, 0.4316, 0.6286));

        observedRegions.add(createObservedRegion(
                chr1.toString(), 9001, 10000, 0.4388, 0.5758));

        return observedRegions;
    }
    
    
    @Test
    public void testNoTumorMode()
    {
        List<ObservedRegion> observedRegions = buildNoTumorObservedRegions();
        
        PurityPloidyFitter fitter;
        fitter = new PurityPloidyFitter(mConfig, mReferenceData, mSampleData, null, mRegionFitCalculator, observedRegions, mSegmentation);
        fitter.run();
        
        PPL_LOGGER.debug(fitter.copyNumberFit());
        PPL_LOGGER.debug(fitter.somaticFit());
        PPL_LOGGER.debug(fitter.finalFit().Method);
        
        assertNotNull(fitter.copyNumberFit());
        
        assertNull(fitter.somaticFit());
        
        String method = fitter.finalFit().Method.toString();
        assertEquals(method, "NO_TUMOR");
    }

    @Test
    public void testCopyNumberMode()
    {
        List<ObservedRegion> observedRegions = buildCopyNumberObservedRegions();

        PurityPloidyFitter fitter;
        fitter = new PurityPloidyFitter(mConfig, mReferenceData, mSampleData, null, mRegionFitCalculator, observedRegions, mSegmentation);
        fitter.run();

        PPL_LOGGER.debug(fitter.copyNumberFit());
        PPL_LOGGER.debug(fitter.somaticFit());
        PPL_LOGGER.debug(fitter.finalFit().Method);

        assertNotNull(fitter.copyNumberFit());

        assertNull(fitter.somaticFit());

        String method = fitter.finalFit().Method.toString();
        assertEquals(method, "NORMAL");
    }
    
    
    
//    @Test
//    public void testSomaticFit()
//    {
//        List<ObservedRegion> observedRegions = buildNoTumorObservedRegions();
//
//        PurityPloidyFitter fitter = new PurityPloidyFitter(mConfig, mReferenceData, mSampleData, null,
//                mRegionFitCalculator, observedRegions, mSegmentation);
//
//        assertTrue(fitter.isValid());
//
//        fitter.run();
//
//        assertNotNull(fitter.finalFit());
//    }

    @Test
    public void testMostDiploidPurity()
    {
        final FittedPurity fp1 = createRandomPurity(0.3, 0.3, 2.3);
        final FittedPurity fp2 = createRandomPurity(0.3, 0.2, 1.9);
        final FittedPurity fp3 = createRandomPurity(0.4, 0.4, 1.8);
        final FittedPurity fp4 = createRandomPurity(0.4, 0.3, 2.05);

        final List<FittedPurity> all = Lists.newArrayList(fp1, fp2, fp3, fp4);
        Collections.shuffle(all);

        final List<FittedPurity> result = BestFit.mostDiploidPerPurity(all);

        assertEquals(2, result.size());
        assertEquals(fp2, result.get(0));
        assertEquals(fp4, result.get(1));
    }

    public static double nextDouble(@NotNull final Random random)
    {
        return Math.round(random.nextDouble() * 10000D) / 10000D;
    }

    @NotNull
    public static ImmutableFittedPurity.Builder createRandomPurityBuilder(@NotNull Random random)
    {
        return ImmutableFittedPurity.builder()
                .purity(nextDouble(random))
                .normFactor(nextDouble(random))
                .score(nextDouble(random))
                .diploidProportion(nextDouble(random))
                .ploidy(nextDouble(random))
                .somaticPenalty(nextDouble(random));
    }

    @NotNull
    private FittedPurity createRandomPurity(double purity, double score, double ploidy)
    {
        Random random = new Random();
        return createRandomPurityBuilder(random).purity(purity).score(score).ploidy(ploidy).build();
    }

}
