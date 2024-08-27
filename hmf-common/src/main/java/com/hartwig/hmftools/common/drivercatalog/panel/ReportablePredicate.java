package com.hartwig.hmftools.common.drivercatalog.panel;

import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.parseAltTranscriptInfo;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

public class ReportablePredicate
{
    private final Map<String, DriverGene> mDriverGeneMap;

    public ReportablePredicate(final DriverCategory type, final List<DriverGene> driverGenes)
    {
        mDriverGeneMap = driverGenes.stream()
                .filter(x -> x.likelihoodType().equals(type) && x.reportSomatic())
                .collect(Collectors.toMap(DriverGene::gene, x -> x));
    }

    public boolean isReportable(final VariantImpact variantImpact, final VariantType type, boolean isHotspot)
    {
        if(isReportable(
            variantImpact.GeneName, type, isHotspot, variantImpact.CanonicalCodingEffect, variantImpact.CanonicalEffect))
        {
            return true;
        }

        if(variantImpact.OtherReportableEffects.isEmpty())
            return false;

        List<AltTranscriptReportableInfo> altTransEffects = parseAltTranscriptInfo(variantImpact.OtherReportableEffects);

        return altTransEffects.stream().anyMatch(x ->
                isReportable(variantImpact.GeneName, type, isHotspot, x.Effect, x.Effects));
    }

    public boolean isReportable(
            final String gene, final VariantType type, boolean isHotspot,
            final CodingEffect codingEffect, final String effects)
    {
        final DriverGene driverGene = mDriverGeneMap.get(gene);

        if(driverGene == null)
            return false;

        if(isHotspot && driverGene.reportSomaticHotspot())
            return true;

        DriverImpact impact = DriverImpact.select(type, codingEffect);

        // splice ranks above missense so if a gene is reportable for missense but not splice, ensure this is handled
        boolean hasMissense = effects.contains(MISSENSE.effect());

        if(hasMissense && driverGene.reportMissenseAndInframe())
            return true;

        switch(impact)
        {
            case NONSENSE:
            case FRAMESHIFT:
                return driverGene.reportNonsenseAndFrameshift();

            case SPLICE:
                return driverGene.reportSplice();

            case MISSENSE:
            case INFRAME:
                return driverGene.reportMissenseAndInframe();

            default:
                return false;
        }
    }
}
