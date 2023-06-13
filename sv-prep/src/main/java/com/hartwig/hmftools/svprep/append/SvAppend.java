package com.hartwig.hmftools.svprep.append;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.append.AppendConstants.BREAKEND_PROXIMITY;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

public class SvAppend
{
    private final AppendConfig mConfig;
    private final Map<String,List<BreakendData>> mChrBreakendMap;

    public SvAppend(final CommandLine cmd)
    {
        mChrBreakendMap = Maps.newHashMap();

        mConfig = new AppendConfig(cmd);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        if(!loadBreakends() || mChrBreakendMap.isEmpty())
            System.exit(1);

        SV_LOGGER.info("SV Append for {} SV breakends", mChrBreakendMap.values().stream().mapToInt(x -> x.size()).sum());
        long startTimeMs = System.currentTimeMillis();

        List<RegionTask> regionTasks = Lists.newArrayList();

        for(Map.Entry<String,List<BreakendData>> entry : mChrBreakendMap.entrySet())
        {
            String chromosome = entry.getKey();

            List<BreakendData> proximateBreakends = Lists.newArrayList();
            int lastPosition = 0;

            for(BreakendData breakend : entry.getValue())
            {
                if(lastPosition > 0 && breakend.Position - lastPosition > BREAKEND_PROXIMITY)
                {
                    regionTasks.add(new RegionTask(mConfig, chromosome, proximateBreakends));
                    proximateBreakends.clear();
                }

                proximateBreakends.add(breakend);
                lastPosition = breakend.Position;
            }

            // add the last
            regionTasks.add(new RegionTask(mConfig, chromosome, proximateBreakends));
        }

        if(mConfig.Threads > 1)
        {
            final List<Callable> callableList = regionTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            regionTasks.forEach(x -> x.call());
        }

        // write appended output VCF
        VcfWriter vcfWriter = new VcfWriter(mConfig);
        vcfWriter.writeBreakends(mChrBreakendMap);

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;

        SV_LOGGER.info("SV Append complete, mins({})", format("%.3f", timeTakenMs / 60000.0));
    }

    private boolean loadBreakends()
    {
        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.InputVcf);

        if(!vcfFileReader.fileValid())
            return false;

        int breakendCount = 0;
        for(VariantContext variant : vcfFileReader.iterator())
        {
            if(!mConfig.SpecificRegions.isEmpty())
            {
                if(mConfig.SpecificRegions.stream().noneMatch(x -> x.containsPosition(variant.getContig(), variant.getStart())))
                    continue;
            }

            String chromosome = mConfig.RefGenVersion.versionedChromosome(variant.getContig());

            List<BreakendData> breakends = mChrBreakendMap.get(chromosome);

            if(breakends == null)
            {
                breakends = Lists.newArrayList();
                mChrBreakendMap.put(chromosome, breakends);
            }

            breakends.add(BreakendData.fromVariant(variant));
            ++breakendCount;
        }

        vcfFileReader.close();

        SV_LOGGER.info("loaded {} SV breakends from {}", breakendCount, mConfig.InputVcf);
        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = AppendConfig.createCmdLineOptions();
        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        SvAppend svAppender = new SvAppend(cmd);
        svAppender.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
