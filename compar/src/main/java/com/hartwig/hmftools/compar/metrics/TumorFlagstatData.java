package com.hartwig.hmftools.compar.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.TUMOR_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_MAPPED_PROPORTION;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class TumorFlagstatData implements ComparableItem
{
    public final BamFlagStats mFlagstat;

    public TumorFlagstatData(final BamFlagStats flagstat)
    {
        mFlagstat = flagstat;
    }

    @Override
    public Category category()
    {
        return TUMOR_FLAGSTAT;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", mFlagstat.mappedProportion()));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final TumorFlagstatData otherData = (TumorFlagstatData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_MAPPED_PROPORTION, mFlagstat.mappedProportion(), otherData.mFlagstat.mappedProportion(), thresholds);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}