package com.hartwig.hmftools.esvee.processor;

import com.hartwig.hmftools.esvee.assembly.AssemblyExtenderCounters;
import com.hartwig.hmftools.esvee.assembly.Counters;
import com.hartwig.hmftools.esvee.assembly.PrimaryAssemblerCounters;
import com.hartwig.hmftools.esvee.util.Counter;

public class OverallCounters extends Counters<OverallCounters>
{
    public final Counter JunctionsProcessed = new Counter("Junctions Processed", false);
    public final Counter PrimaryAssemblyTime = new Counter("Primary Assembly Time", true);
    public final Counter InterJunctionDeduplicationTime = new Counter("Deduplication Time", true);
    public final Counter ExtensionTime = new Counter("Extension Time", true);
    public final Counter PrimaryPhasingTime = new Counter("Primary Phasing Time", true);
    public final Counter PhasedAssemblyMergingTime = new Counter("Phased Assembly Merge Time", true);
    public final Counter SecondaryPhasingTime = new Counter("Secondary Phasing Time", true);
    public final Counter MergeSecondaryTime = new Counter("Merge Secondary Time", true);
    public final Counter AlignmentTime = new Counter("Alignment Time", true);
    public final Counter HomologyTime = new Counter("Homology Time", true);
    public final Counter VariantCallingTime = new Counter("Variant Calling Time", true);
    public final Counter SupportScanTime = new Counter("Rescan Support Time", true);

    public final Counter ExtraScannedSupport = new Counter("Extra Support", false);

    public final com.hartwig.hmftools.esvee.assembly.PrimaryAssemblerCounters PrimaryAssemblerCounters = new PrimaryAssemblerCounters();
    public final com.hartwig.hmftools.esvee.assembly.AssemblyExtenderCounters AssemblyExtenderCounters = new AssemblyExtenderCounters();
    public final VariantDeduplicationCounters VariantDeduplicationCounters = new VariantDeduplicationCounters();
}