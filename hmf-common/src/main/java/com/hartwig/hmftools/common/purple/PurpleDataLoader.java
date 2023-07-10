package com.hartwig.hmftools.common.purple;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.AllTranscriptSomaticVariant;
import com.hartwig.hmftools.common.variant.AllTranscriptSomaticVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public final class PurpleDataLoader
{
    private PurpleDataLoader()
    {
    }

    @NotNull
    public static PurpleData load(final String tumorSample, @Nullable final String referenceSample, @Nullable final String rnaSample,
            final String purpleDir) throws IOException
    {
        String qcFile = PurpleQCFile.generateFilename(purpleDir, tumorSample);
        String purityTsv = PurityContextFile.generateFilenameForReading(purpleDir, tumorSample);
        String somaticDriverCatalogTsv = DriverCatalogFile.generateSomaticFilename(purpleDir, tumorSample);
        String somaticVariantVcf = resolveVcfPath(PurpleCommon.purpleSomaticVcfFile(purpleDir, tumorSample));
        String germlineDriverCatalogTsv = DriverCatalogFile.generateGermlineFilename(purpleDir, tumorSample);
        String germlineVariantVcf = resolveVcfPath(PurpleCommon.purpleGermlineVcfFile(purpleDir, tumorSample));
        String somaticStructuralVariantVcf = resolveVcfPath(PurpleCommon.purpleSomaticSvFile(purpleDir, tumorSample));
        String germlineStructuralVariantVcf = resolveVcfPath(PurpleCommon.purpleGermlineSvFile(purpleDir, tumorSample));
        String copyNumberTsv = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, tumorSample);
        String geneCopyNumberTsv = GeneCopyNumberFile.generateFilename(purpleDir, tumorSample);
        String germlineDeletionTsv = GermlineDeletion.generateFilename(purpleDir, tumorSample);

        return load(tumorSample,
                referenceSample,
                rnaSample,
                qcFile,
                purityTsv,
                somaticDriverCatalogTsv,
                somaticVariantVcf,
                germlineDriverCatalogTsv,
                germlineVariantVcf,
                somaticStructuralVariantVcf,
                germlineStructuralVariantVcf,
                copyNumberTsv,
                geneCopyNumberTsv,
                germlineDeletionTsv);
    }

    private static String resolveVcfPath(final String vcfPath)
    {
        if (!new File(vcfPath).exists() && vcfPath.endsWith(".gz"))
        {
            String unzippedVcfPath = vcfPath.substring(0, vcfPath.length() - 3);
            if (new File(unzippedVcfPath).exists())
            {
                return unzippedVcfPath;
            }
        }
        return vcfPath;
    }

    @NotNull
    private static PurpleData load(@NotNull String tumorSample, @Nullable String referenceSample, @Nullable String rnaSample,
            @NotNull String qcFile, @NotNull String purityTsv, @NotNull String somaticDriverCatalogTsv, @NotNull String somaticVariantVcf,
            @NotNull String germlineDriverCatalogTsv, @NotNull String germlineVariantVcf, @NotNull String somaticStructuralVariantVcf,
            @NotNull String germlineStructuralVariantVcf, @NotNull String copyNumberTsv, @NotNull String geneCopyNumberTsv,
            @NotNull String germlineDeletionTsv) throws IOException
    {
        PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityTsv);

        List<DriverCatalog> somaticDrivers = DriverCatalogFile.read(somaticDriverCatalogTsv);

        List<AllTranscriptSomaticVariant> allSomaticVariants =
                AllTranscriptSomaticVariantFactory.passOnlyInstance().fromVCFFile(tumorSample, referenceSample, rnaSample, somaticVariantVcf);
        List<AllTranscriptSomaticVariant> reportableSomaticVariants = selectReportedVariants(allSomaticVariants);

        List<StructuralVariant> allSomaticStructuralVariants =
                StructuralVariantFileLoader.fromFile(somaticStructuralVariantVcf, new PassingVariantFilter());

        List<GeneCopyNumber> allSomaticGeneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberTsv);

        List<DriverCatalog> germlineDrivers = null;
        List<StructuralVariant> allGermlineStructuralVariants = null;
        List<AllTranscriptSomaticVariant> allGermlineVariants = null;
        List<AllTranscriptSomaticVariant> reportableGermlineVariants = null;
        List<GermlineDeletion> allGermlineDeletions = null;
        List<GermlineDeletion> reportableGermlineDeletions = null;
        if (referenceSample != null)
        {
            germlineDrivers = DriverCatalogFile.read(germlineDriverCatalogTsv);
            allGermlineStructuralVariants = StructuralVariantFileLoader.fromFile(germlineStructuralVariantVcf, new PassingVariantFilter());

            allGermlineVariants = new AllTranscriptSomaticVariantFactory().fromVCFFile(tumorSample, referenceSample, rnaSample, germlineVariantVcf);
            reportableGermlineVariants = selectReportedVariants(allGermlineVariants);

            allGermlineDeletions = selectPassDeletions(GermlineDeletion.read(germlineDeletionTsv));
            reportableGermlineDeletions = selectReportedDeletions(allGermlineDeletions);
        }

        return ImmutablePurpleData.builder()
                .purityContext(purityContext)
                .somaticDrivers(somaticDrivers)
                .germlineDrivers(germlineDrivers)
                .allSomaticVariants(allSomaticVariants)
                .reportableSomaticVariants(reportableSomaticVariants)
                .allGermlineVariants(allGermlineVariants)
                .reportableGermlineVariants(reportableGermlineVariants)
                .allSomaticStructuralVariants(allSomaticStructuralVariants)
                .allGermlineStructuralVariants(allGermlineStructuralVariants)
                .allSomaticCopyNumbers(PurpleCopyNumberFile.read(copyNumberTsv))
                .allSomaticGeneCopyNumbers(allSomaticGeneCopyNumbers)
                .allGermlineDeletions(allGermlineDeletions)
                .reportableGermlineDeletions(reportableGermlineDeletions)
                .build();
    }

    @NotNull
    private static List<AllTranscriptSomaticVariant> selectReportedVariants(@NotNull List<AllTranscriptSomaticVariant> allVariants)
    {
        List<AllTranscriptSomaticVariant> reported = Lists.newArrayList();
        for (AllTranscriptSomaticVariant variant : allVariants)
        {
            if (variant.reported())
            {
                reported.add(variant);
            }
        }
        return reported;
    }

    @NotNull
    private static List<GermlineDeletion> selectPassDeletions(@NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        List<GermlineDeletion> pass = Lists.newArrayList();
        for (GermlineDeletion deletion : allGermlineDeletions)
        {
            if (deletion.Filter.equals("PASS"))
            {
                pass.add(deletion);
            }
        }

        return pass;
    }

    @NotNull
    private static List<GermlineDeletion> selectReportedDeletions(@NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        List<GermlineDeletion> reported = Lists.newArrayList();
        for (GermlineDeletion deletion : allGermlineDeletions)
        {
            if (deletion.Reported)
            {
                reported.add(deletion);
            }
        }
        return reported;
    }
}
