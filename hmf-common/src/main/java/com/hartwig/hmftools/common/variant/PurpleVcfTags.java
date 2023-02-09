package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class PurpleVcfTags
{
    public static final String PURPLE_JUNCTION_COPY_NUMBER = "PURPLE_JCN";
    public static final String PURPLE_JUNCTION_COPY_NUMBER_DESC = "Purity adjusted copy number of variant junction";

    public static final String PURPLE_PLOIDY_INFO = "PURPLE_PLOIDY";

    public static final String PURPLE_CN_CHANGE = "PURPLE_CN_CHANGE";
    public static final String PURPLE_CN_CHANGE_DESC = "Purity adjusted change in copy number at each breakend";

    public static final String PURPLE_AF = "PURPLE_AF";
    public static final String PURPLE_AF_DESC = "Purity adjusted variant allelic frequency";

    public static final String PURPLE_CN = "PURPLE_CN";
    public static final String PURPLE_CN_DESC = "Purity adjusted copy number surrounding variant location";

    public static final String PURPLE_BIALLELIC_FLAG = "BIALLELIC";
    public static final String PURPLE_BIALLELIC_DESC = "Variant is biallelic";

    public static final String PURPLE_VARIANT_CN = "PURPLE_VCN";
    private static final String PURPLE_VARIANT_CN_DESC = "Purity adjusted variant copy number";

    public static final String PURPLE_GERMLINE_INFO = "PURPLE_GERMLINE";
    private static final String PURPLE_GERMLINE_DESC = "Germline classification surrounding variant location";

    public static final String PURPLE_MINOR_ALLELE_CN_INFO = "PURPLE_MACN";
    private static final String PURPLE_MINOR_ALLELE_PLOIDY_DESC = "Purity adjusted minor allele ploidy surrounding variant location";

    public static VCFHeader addGermlineHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template)
    {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN, 1, VCFHeaderLineType.Float, PURPLE_VARIANT_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));

        return template;
    }

    @NotNull
    public static VCFHeader addSomaticHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template)
    {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN, 1, VCFHeaderLineType.Float, PURPLE_VARIANT_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_GERMLINE_INFO, 1, VCFHeaderLineType.String, PURPLE_GERMLINE_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));

        return template;
    }
}