package com.hartwig.hmftools.common.purple;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.variant.AllTranscriptSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleData
{
    @NotNull
    PurityContext purityContext();

    @NotNull
    List<DriverCatalog> somaticDrivers();

    @Nullable
    List<DriverCatalog> germlineDrivers();

    @NotNull
    List<AllTranscriptSomaticVariant> allSomaticVariants();

    @NotNull
    List<AllTranscriptSomaticVariant> reportableSomaticVariants();

    @Nullable
    List<AllTranscriptSomaticVariant> allGermlineVariants();

    @Nullable
    List<AllTranscriptSomaticVariant> reportableGermlineVariants();

    @NotNull
    List<StructuralVariant> allSomaticStructuralVariants();

    @Nullable
    List<StructuralVariant> allGermlineStructuralVariants();

    @NotNull
    List<PurpleCopyNumber> allSomaticCopyNumbers();

    @NotNull
    List<GeneCopyNumber> allSomaticGeneCopyNumbers();

    @Nullable
    List<GermlineDeletion> allGermlineDeletions();

    @Nullable
    List<GermlineDeletion> reportableGermlineDeletions();

}
