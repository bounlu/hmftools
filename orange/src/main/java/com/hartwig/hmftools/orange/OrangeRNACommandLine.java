package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.orange.OrangeCommandLine.ORANGE_JSON;

import java.io.IOException;
import java.util.Optional;

import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRNAConfig;
import com.hartwig.hmftools.datamodel.orange.OrangeRNAConfig;
import com.hartwig.hmftools.orange.util.Config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public interface OrangeRNACommandLine {

    Logger LOGGER = LogManager.getLogger(OrangeRNACommandLine.class);

    String RNA_SAMPLE_ID = "rna_sample_id";

    String ISOFOX_GENE_DISTRIBUTION_CSV = "isofox_gene_distribution_csv";
    String ISOFOX_ALT_SJ_COHORT_CSV = "isofox_alt_sj_cohort_csv";

    String ISOFOX_SUMMARY_CSV = "isofox_summary_csv";
    String ISOFOX_GENE_DATA_CSV = "isofox_gene_data_csv";
    String ISOFOX_FUSION_CSV = "isofox_fusion_csv";
    String ISOFOX_ALT_SPLICE_JUNCTION_CSV = "isofox_alt_splice_junction_csv";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        options.addOption(ORANGE_JSON, true, "(Optional) Location of an existing orange json");

        options.addOption(RNA_SAMPLE_ID, true, "(Optional) The RNA sample of the tumor sample for which ORANGE will run.");

        options.addOption(ISOFOX_GENE_DISTRIBUTION_CSV, true, "(Optional) Path to isofox gene distribution CSV.");
        options.addOption(ISOFOX_ALT_SJ_COHORT_CSV, true, "(Optional) Path to isofox alt SJ cohort CSV.");

        options.addOption(ISOFOX_SUMMARY_CSV, true, "(Optional) Path towards the ISOFOX summary data.");
        options.addOption(ISOFOX_GENE_DATA_CSV, true, "(Optional) Path towards the ISOFOX gene data.");
        options.addOption(ISOFOX_FUSION_CSV, true, "(Optional) Path towards the ISOFOX fusion data.");
        options.addOption(ISOFOX_ALT_SPLICE_JUNCTION_CSV, true, "(Optional) Path towards the ISOFOX alt splice junction data.");

        return options;
    }

    @NotNull
    String rnaSampleId();

    @NotNull
    String isofoxGeneDistributionCsv();

    @NotNull
    String isofoxAltSjCohortCsv();

    @NotNull
    String isofoxSummaryCsv();

    @NotNull
    String isofoxGeneDataCsv();

    @NotNull
    String isofoxFusionCsv();

    @NotNull
    String isofoxAltSpliceJunctionCsv();

    @NotNull
    static OrangeRNAConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        if (!hasCompleteRnaConfig(cmd,
                RNA_SAMPLE_ID,
                ISOFOX_GENE_DISTRIBUTION_CSV,
                ISOFOX_ALT_SJ_COHORT_CSV,
                ISOFOX_SUMMARY_CSV,
                ISOFOX_GENE_DATA_CSV,
                ISOFOX_FUSION_CSV,
                ISOFOX_ALT_SPLICE_JUNCTION_CSV)) {
            LOGGER.debug("No proper RNA input fed to ORANGE");
            return null;
        }

        if (cmd.hasOption(ORANGE_JSON)) {
            try {
                return Optional.ofNullable(OrangeJson.getInstance().read(cmd.getOptionValue(ORANGE_JSON)).config().rnaConfig())
                        .orElseThrow(() -> new IllegalStateException("Cannot rerun ORANGE rna from a JSON without RNA configuration"));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        String rnaSampleId = Config.nonOptionalValue(cmd, RNA_SAMPLE_ID);
        LOGGER.debug("RNA sample configured as {}", rnaSampleId);

        return ImmutableOrangeRNAConfig.builder()
                .rnaSampleId(rnaSampleId)
                .isofoxGeneDistributionCsv(Config.nonOptionalFile(cmd, ISOFOX_GENE_DISTRIBUTION_CSV))
                .isofoxAltSjCohortCsv(Config.nonOptionalFile(cmd, ISOFOX_ALT_SJ_COHORT_CSV))
                .isofoxSummaryCsv(Config.nonOptionalFile(cmd, ISOFOX_SUMMARY_CSV))
                .isofoxGeneDataCsv(Config.nonOptionalFile(cmd, ISOFOX_GENE_DATA_CSV))
                .isofoxFusionCsv(Config.nonOptionalFile(cmd, ISOFOX_FUSION_CSV))
                .isofoxAltSpliceJunctionCsv(Config.nonOptionalFile(cmd, ISOFOX_ALT_SPLICE_JUNCTION_CSV))
                .build();
    }

    private static boolean hasCompleteRnaConfig(@NotNull CommandLine cmd, @NotNull String... options) {
        boolean hasAllOptions = true;
        boolean hasAnyOption = false;
        for (String option : options) {
            boolean hasOption = cmd.hasOption(option);
            hasAllOptions = hasAllOptions && hasOption;
            hasAnyOption = hasAnyOption || hasOption;
        }

        if (hasAnyOption && !hasAllOptions) {
            LOGGER.warn("RNA config has been provided but incompletely");
        }

        return hasAllOptions;
    }
}
