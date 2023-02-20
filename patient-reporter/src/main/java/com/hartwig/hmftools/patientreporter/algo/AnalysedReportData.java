package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.patientreporter.ReportData;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.remarks.SpecialRemarkModel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AnalysedReportData implements ReportData {

    @NotNull
    public abstract GermlineReportingModel germlineReportingModel();

    @NotNull
    public abstract SpecialRemarkModel specialRemarkModel();

    @NotNull
    public abstract KnownFusionCache knownFusionCache();
}