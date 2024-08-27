package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;

import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_READ_MAX_MISMATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.buildIndelFrequencies;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensionReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findMaxFrequencyIndelReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.hasIndelJunctionReads;
import static com.hartwig.hmftools.esvee.assembly.SequenceCompare.calcReadSequenceMismatches;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.recordSoftClipsAtJunction;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.assembly.types.ReadAssemblyIndices.getJunctionReadExtensionIndices;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.types.ReadAssemblyIndices;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadFilters;

public class JunctionAssembler
{
    private final Junction mJunction;
    private final List<Read> mNonJunctionReads;

    public JunctionAssembler(final Junction junction)
    {
        mJunction = junction;
        mNonJunctionReads = Lists.newArrayList();
    }

    public List<Read> nonJunctionReads() { return mNonJunctionReads; }

    public List<JunctionAssembly> processJunction(final List<Read> rawReads)
    {
        // find prominent reads to establish the extension sequence, taking any read meeting min soft-clip lengths
        // and repetitive indels

        List<Read> junctionReads = Lists.newArrayList();
        List<Read> extensionReads = Lists.newArrayList();

        if(!mJunction.indelBased() && hasIndelJunctionReads(mJunction, rawReads))
        {
            // fall-back in case Prep didn't set this state or junctions are loaded from config
            mJunction.markAsIndel();
        }

        boolean hasMinLengthSoftClipRead = false;

        if(mJunction.indelBased())
        {
            findIndelExtensionReads(mJunction, rawReads, extensionReads, junctionReads, mNonJunctionReads);
            hasMinLengthSoftClipRead = !extensionReads.isEmpty();
        }
        else
        {
            Map<Integer,List<Read>> indelLengthReads = Maps.newHashMap();

            // the only difference for indel-based junctions is that only the long indels are used to build the consensus extension

            for(Read read : rawReads)
            {
                if(!ReadFilters.recordSoftClipsAndCrossesJunction(read, mJunction))
                {
                    mNonJunctionReads.add(read);
                    continue;
                }

                if(recordSoftClipsAtJunction(read, mJunction))
                {
                    int softClipJunctionExtension = readJunctionExtensionLength(read, mJunction);

                    hasMinLengthSoftClipRead |= softClipJunctionExtension >= ASSEMBLY_MIN_SOFT_CLIP_LENGTH
                            || (read.hasLineTail() && softClipJunctionExtension >= LINE_MIN_EXTENSION_LENGTH);

                    if(softClipJunctionExtension >= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH)
                    {
                        extensionReads.add(read);
                    }
                }

                if((mJunction.isForward() && read.indelImpliedAlignmentEnd() > 0)
                || (mJunction.isReverse() && read.indelImpliedAlignmentStart() > 0))
                {
                    buildIndelFrequencies(indelLengthReads, read);
                }

                junctionReads.add(read);
            }

            List<Read> dominantIndelReads = findMaxFrequencyIndelReads(indelLengthReads);

            extensionReads.addAll(dominantIndelReads);
        }

        if(!hasMinLengthSoftClipRead || extensionReads.size() < ASSEMBLY_MIN_READ_SUPPORT)
            return Collections.emptyList();

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);

        int reqExtensionLength = extensionSeqBuilder.hasLineSequence() ? LINE_MIN_EXTENSION_LENGTH : ASSEMBLY_MIN_SOFT_CLIP_LENGTH;

        if(!extensionSeqBuilder.isValid() || extensionSeqBuilder.extensionLength() < reqExtensionLength)
            return Collections.emptyList();

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        if(assemblySupport.size() < ASSEMBLY_MIN_READ_SUPPORT)
            return Collections.emptyList();

        JunctionAssembly firstAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        if(extensionSeqBuilder.hasLineSequence())
            firstAssembly.markLineSequence();

        List<JunctionAssembly> assemblies = Lists.newArrayList(firstAssembly);

        // test for a second well-supported, alternative assembly at the same junction
        JunctionAssembly secondAssembly = checkSecondAssembly(extensionSeqBuilder.mismatchReads(), firstAssembly);

        if(secondAssembly != null)
            assemblies.add(secondAssembly);

        for(JunctionAssembly assembly : assemblies)
        {
            int mismatchReadCount = 0;

            // test other junction-spanning reads against this new assembly
            for(Read read : junctionReads)
            {
                if(assembly.support().stream().anyMatch(x -> x.cachedRead() == read)) // skip those already added
                    continue;

                if(!canAddJunctionRead(assembly, read))
                    ++mismatchReadCount;
            }

            assembly.addMismatchReadCount(mismatchReadCount);

            RefBaseSeqBuilder refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);
            assembly.setRefBases(refBaseSeqBuilder);

            assembly.buildRepeatInfo();
        }

        return assemblies;
    }

    private JunctionAssembly checkSecondAssembly(final List<Read> extensionReads, final JunctionAssembly firstAssembly)
    {
        if(extensionReads.isEmpty())
            return null;

        if(firstAssembly.hasLineSequence())
            return null;

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);

        if(!extensionSeqBuilder.isValid())
            return null;

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        int secondSupport = assemblySupport.size();
        double secondSupportPerc = secondSupport / (double)firstAssembly.supportCount();

        if(secondSupport < ASSEMBLY_SPLIT_MIN_READ_SUPPORT || secondSupportPerc < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC)
            return null;

        JunctionAssembly newAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        if(newAssembly.extensionLength() < ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
            return null;

        // perform a final sequence comparison check with more liberal comparison tests
        boolean closeMatch = SequenceCompare.matchedAssemblySequences(firstAssembly, newAssembly);
        return !closeMatch ? newAssembly : null;
    }

    private boolean canAddJunctionRead(final JunctionAssembly assembly, final Read read)
    {
        int readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);

        if(readJunctionIndex == INVALID_INDEX)
            return false;

        ReadAssemblyIndices readAssemblyIndices = getJunctionReadExtensionIndices(
                assembly.junction(), assembly.junctionIndex(), read, readJunctionIndex);

        int assemblyIndexStart = readAssemblyIndices.AssemblyIndexStart;
        int readIndexStart = readAssemblyIndices.ReadIndexStart;
        int readIndexEnd = readAssemblyIndices.ReadIndexEnd;

        if(assemblyIndexStart < 0)
        {
            // allow for indel-adjusted reads
            if(read.indelImpliedAlignmentStart() != mJunction.Position)
                return false;

            readIndexStart -= assemblyIndexStart;
            assemblyIndexStart = 0;
        }

        // first attempt a straight string match for simplicity
        int matchLength = readIndexEnd - readIndexStart + 1;

        if(matchLength < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
            return false;

        int highQualMatchCount = 0;
        int mismatchCount = 0;
        int checkedBaseCount = 0;

        final byte[] assemblyBases = assembly.bases();
        final byte[] assemblyBaseQuals = assembly.baseQuals();

        int assemblyIndex = assemblyIndexStart;
        int assemblyBaseLength = assembly.baseLength();

        for(int i = readIndexStart; i <= readIndexEnd; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= assemblyBaseLength)
                break;

            byte qual = read.getBaseQuality()[i];
            ++checkedBaseCount;

            if(basesMatch(read.getBases()[i], assemblyBases[assemblyIndex], qual, assemblyBaseQuals[assemblyIndex]))
            {
                if(aboveMinQual(qual) && assemblyIndex != assembly.junctionIndex())
                    ++highQualMatchCount;
            }
            else
            {
                ++mismatchCount;

                if(mismatchCount > PRIMARY_ASSEMBLY_READ_MAX_MISMATCH)
                    break;
            }
        }

        int permittedMismatches = mismatchesPerComparisonLength(checkedBaseCount);

        if(mismatchCount > permittedMismatches)
        {
            checkedBaseCount = readIndexEnd - readIndexStart + 1;

            if(assemblyIndex < 0)
                checkedBaseCount = max(checkedBaseCount + assemblyIndex, 0);

            permittedMismatches = mismatchesPerComparisonLength(checkedBaseCount);

            // test again taking repeats into consideration
            mismatchCount = calcReadSequenceMismatches(
                    mJunction.isForward(), assemblyBases, assemblyBaseQuals, assembly.repeatInfo(), read, readJunctionIndex, permittedMismatches);
        }

        if(mismatchCount > permittedMismatches || highQualMatchCount < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
            return false;

        assembly.addSupport(read, JUNCTION, readJunctionIndex, highQualMatchCount, mismatchCount);
        return true;
    }
}
