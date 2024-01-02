package com.hartwig.hmftools.sage.evidence;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.candidate_.Candidate_;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.quality.QualityCalculator;

public class ReadContextCounterFactory
{
    private static final Set<VariantTier> HIGH_COVERAGE = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    private final SageConfig mConfig;

    public ReadContextCounterFactory(final SageConfig config)
    {
        mConfig = config;
    }

    public List<ReadContextCounter> create(
            final List<Candidate_> candidates, final SageConfig config, final QualityCalculator qualityCalculator, final String sampleId)
    {
        List<ReadContextCounter> readCounters = Lists.newArrayListWithExpectedSize(candidates.size());

        int readId = 0;

        for(Candidate_ candidate : candidates)
        {
            readCounters.add(new ReadContextCounter(
                    readId++,
                    candidate.variant(), candidate.readContext(), candidate.tier(),
                    maxCoverage(candidate), candidate.minNumberOfEvents(),
                    config, qualityCalculator, sampleId));
        }

        return readCounters;
    }

    private int maxCoverage(final Candidate_ candidate)
    {
        return HIGH_COVERAGE.contains(candidate.tier()) ? mConfig.MaxReadDepthPanel : mConfig.MaxReadDepth;
    }
}
