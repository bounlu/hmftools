package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogMapTest;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.test.SomaticVariantTestBuilderFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canResolveReportableSomaticVariants() {
        String gene1 = "gene1";
        String gene2 = "gene2";
        SomaticVariant variant1 =
                SomaticVariantTestBuilderFactory.create().reported(true).gene(gene1).canonicalTranscript("transcript1").build();
        SomaticVariant variant2 =
                SomaticVariantTestBuilderFactory.create().reported(false).gene(gene2).canonicalTranscript("transcript2").build();

        double likelihood = 0.6;
        DriverCatalog driverGene1 = DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript1");

        List<ReportableVariant> reportable = ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant1, variant2),
                Lists.newArrayList(driverGene1));

        assertEquals(1, reportable.size());
        assertEquals(likelihood, reportable.get(0).driverLikelihood(), EPSILON);
    }

    @Test
    public void canResolveGermlineVariantsWithMultipleDrivers() {
        String gene = "gene";
        SomaticVariant variant =
                SomaticVariantTestBuilderFactory.create().reported(true).gene(gene).canonicalTranscript("transcript1").build();

        DriverCatalog driver1 = DriverCatalogMapTest.createCanonicalGermlineMutationEntryForGene(gene, 0.6, "transcript1");
        DriverCatalog driver2 =
                ImmutableDriverCatalog.builder().from(driver1).driver(DriverType.GERMLINE_DELETION).driverLikelihood(1D).build();

        List<ReportableVariant> reportable =
                ReportableVariantFactory.toReportableGermlineVariants(Lists.newArrayList(variant), Lists.newArrayList(driver1, driver2));

        assertEquals(0.6, reportable.get(0).driverLikelihood(), EPSILON);
    }

    @Test
    public void canResolveReportableFromNonCanonicalDrivers() {
        String gene = "gene";
        SomaticVariant variant = SomaticVariantTestBuilderFactory.create()
                .reported(true)
                .gene(gene)
                .canonicalTranscript("transcript1")
                .otherReportedEffects("ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT")
                .build();

        double likelihood = 0.6;
        DriverCatalog driverNonCanonical = ImmutableDriverCatalog.builder()
                .from(DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene, likelihood, "ENST00000579755"))
                .isCanonical(false)
                .build();

        List<ReportableVariant> reportable =
                ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant), Lists.newArrayList(driverNonCanonical));

        assertEquals(1, reportable.size());
        assertEquals(likelihood, reportable.get(0).driverLikelihood(), EPSILON);
        //assertEquals("c.246_247delCG", reportable.get(0).canonicalHgvsCodingImpact());
       // assertEquals("p.Gly83fs", reportable.get(0).canonicalHgvsProteinImpact());

        double likelihoodCanonical = 0.6;
        SomaticVariant variant2 = SomaticVariantTestBuilderFactory.create()
                .reported(true)
                .gene(gene)
                .canonicalTranscript("transcript2")
                .otherReportedEffects("ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT")
                .build();
        DriverCatalog driverCanonical = ImmutableDriverCatalog.builder()
                .from(DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene, likelihoodCanonical, "transcript2"))
                .isCanonical(true)
                .build();
        List<ReportableVariant> reportable2 = ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant2),
                Lists.newArrayList(driverNonCanonical, driverCanonical));

        assertEquals(2, reportable2.size());

//        ReportableVariant reportableVariant1 = extractVariant(true, reportable2);
//        assertEquals(likelihoodCanonical, reportableVariant1.driverLikelihood(), EPSILON);
//        assertEquals(Strings.EMPTY, reportableVariant1.canonicalHgvsCodingImpact());
//        assertEquals(Strings.EMPTY, reportableVariant1.canonicalHgvsProteinImpact());
//
//        ReportableVariant reportableVariant2 = extractVariant(false, reportable2);
//        assertEquals(likelihoodCanonical, reportableVariant2.driverLikelihood(), EPSILON);
//        assertEquals(Strings.EMPTY, reportableVariant2.canonicalHgvsCodingImpact());
//        assertEquals(Strings.EMPTY, reportableVariant2.canonicalHgvsProteinImpact());

    }

    @NotNull
    public static ReportableVariant extractVariant(boolean canonical, @NotNull List<ReportableVariant> reportableVariants) {
        for (ReportableVariant reportableVariant: reportableVariants) {
            if (reportableVariant.isCanonical() && canonical) {
                return reportableVariant;
            }
        }
        throw new IllegalStateException("Could not find reportbale variant : " + canonical);
    }
}