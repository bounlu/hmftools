package com.hartwig.hmftools.pave.impact;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Codons.CODON_LENGTH;
import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.pave.impact.ImpactClassifier.checkStopStartCodons;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.impact.ProteinUtils.trimAminoAcids;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.pave.VariantData;

public class PhasedVariantClassifier
{
    private final List<PhasedVariants> mPhasedVariants;

    public PhasedVariantClassifier()
    {
        mPhasedVariants = Lists.newArrayList();
    }

    public void checkAddVariant(final VariantData variant)
    {
        if(!variant.hasLocalPhaseSet())
            return;

        PhasedVariants phasedVariants = mPhasedVariants.stream()
                .filter(x -> x.LocalPhaseId == variant.localPhaseSet()).findFirst().orElse(null);

        if(phasedVariants == null)
        {
            phasedVariants = new PhasedVariants(variant.localPhaseSet());
            mPhasedVariants.add(phasedVariants);
        }

        phasedVariants.addVariant(variant);
    }

    public boolean hasCompleteVariants(int currentLocalPhaseSet)
    {
        return !mPhasedVariants.isEmpty() && mPhasedVariants.get(0).LocalPhaseId != currentLocalPhaseSet;
    }

    public List<PhasedVariants> popCompletePhasedVariants(int currentLocalPhaseSet)
    {
        List<PhasedVariants> completeVariants = Lists.newArrayList();

        while(!mPhasedVariants.isEmpty())
        {
            if(mPhasedVariants.get(0).LocalPhaseId != currentLocalPhaseSet)
            {
                completeVariants.add(mPhasedVariants.get(0));
                mPhasedVariants.remove(0);
            }
            else
            {
                break;
            }
        }

        return completeVariants;
    }

    public void clear() { mPhasedVariants.clear(); }

    public void reclassifyPhasedVariants(final PhasedVariants phasedVariants, final RefGenomeInterface refGenome)
    {
        if(phasedVariants.variants().size() < 2)
            return;

        // ignore groups without INDELS
        if(phasedVariants.variants().stream().noneMatch(x -> x.isIndel()))
            return;

        // PV_LOGGER.trace("processing {} phased variant variants with localPhaseId({})", variants.Variants.size(), variants.LocalPhaseId);

        // assume that variants are ordered  by ascending position
        VariantData variant = phasedVariants.variants().get(0);

        if(variant.getImpacts().isEmpty())
            return;

        // process impacts by matching gene and transcript
        for(Map.Entry<String,List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
        {
            String gene = entry.getKey();

            for(VariantTransImpact transImpact : entry.getValue())
            {
                List<VariantTransImpact> transImpacts = Lists.newArrayList(transImpact);
                List<VariantData> variants = Lists.newArrayList(variant);

                for(int i = 1; i < phasedVariants.variants().size(); ++i)
                {
                    VariantData nextVariant = phasedVariants.variants().get(i);

                    List<VariantTransImpact> otherImpacts = nextVariant.getImpacts().get(gene);

                    if(otherImpacts != null)
                    {
                        VariantTransImpact nextImpact = otherImpacts.stream()
                                .filter(x -> x.TransData.TransName.equals(transImpact.TransData.TransName))
                                .findFirst().orElse(null);

                        if(nextImpact != null)
                        {
                            transImpacts.add(nextImpact);
                            variants.add(nextVariant);
                        }
                    }
                }

                reclassifyImpacts(phasedVariants.LocalPhaseId, variants, transImpacts, refGenome);
            }
        }
    }

    private void reclassifyImpacts(
            int localPhaseSet, final List<VariantData> variants, final List<VariantTransImpact> transImpacts, final RefGenomeInterface refGenome)
    {
        // ignore if not 2 or more with coding impacts
        if(transImpacts.stream().filter(x -> x.hasCodingBases()).count() < 2)
            return;

        if(transImpacts.stream().anyMatch(x -> x.hasEffect(STOP_LOST) || x.hasEffect(START_LOST)))
            return;

        // ignore if not phased anyway
        int indelBaseTotal = 0;
        int frameshiftCount = 0;
        int minIndelIndex = transImpacts.size();
        int maxIndelIndex = 0;

        for(int i = 0; i < transImpacts.size(); ++i)
        {
            VariantTransImpact transImpact = transImpacts.get(i);
            VariantData variant = variants.get(i);

            if(!transImpact.hasCodingBases() || !variant.isIndel())
                continue;

            indelBaseTotal += variant.isInsert() ? variant.baseDiff() : -transImpact.codingContext().DeletedCodingBases;

            if(transImpact.codingContext().IsFrameShift)
                ++frameshiftCount;

            minIndelIndex = min(minIndelIndex, i);
            maxIndelIndex = max(maxIndelIndex, i);
        }

        // ignore any SNV or MNV which is not in between INDELs and doesn't overlap another variant
        List<VariantData> ignoredVariants = Lists.newArrayList();
        Map<VariantTransImpact,ImpactedRefCodingData> refCodingImpacts = Maps.newHashMap();

        // at the same time look for an overlap between an inframe INDEL and any SNV/MNV
        boolean hasInframeSnvMnvOverlap = false;

        for(int i = 0; i < transImpacts.size(); ++i)
        {
            VariantTransImpact transImpact = transImpacts.get(i);
            VariantData variant = variants.get(i);

            if(i > 0 && variants.get(i - 1).Position >= variant.Position)
                return; // don't handle overlapping or out-of-order variants

            if(!transImpact.hasCodingBases() || variant.isIndel())
                continue;

            VariantTransImpact prevTransImpact = i > 0 ? transImpacts.get(i - 1) : null;

            boolean overlapsOnStart = false;

            ImpactedRefCodingData transRefCodingData = getImpactedRefCodingData(refCodingImpacts, transImpact);

            if(prevTransImpact != null && prevTransImpact.hasCodingBases())
            {
                ImpactedRefCodingData prevRefCodingData = getImpactedRefCodingData(refCodingImpacts, prevTransImpact);
                overlapsOnStart = prevRefCodingData.PosEnd >= transRefCodingData.PosStart;
            }

            boolean overlapsOnEnd = false;
            VariantTransImpact nextTransImpact = null;

            if(!overlapsOnStart)
            {
                nextTransImpact = i < transImpacts.size() - 1 ? transImpacts.get(i + 1) : null;
                if(nextTransImpact != null && nextTransImpact.hasCodingBases())
                {
                    ImpactedRefCodingData nextRefCodingData = getImpactedRefCodingData(refCodingImpacts, nextTransImpact);
                    overlapsOnEnd = transRefCodingData.PosEnd >= nextRefCodingData.PosStart;
                }
            }

            if(i < minIndelIndex || i > maxIndelIndex)
            {
                if(!overlapsOnStart && !overlapsOnEnd)
                {
                    // an SNV/MNV outside the INDEL range and not overlapping
                    ignoredVariants.add(variant);
                    continue;
                }
            }

            // otherwise check if the overlap is with an inframe INDEL
            if(overlapsOnStart && variants.get(i - 1).isIndel() && !prevTransImpact.codingContext().IsFrameShift)
                hasInframeSnvMnvOverlap = true;

            if(overlapsOnEnd && variants.get(i + 1).isIndel() && !nextTransImpact.codingContext().IsFrameShift)
                hasInframeSnvMnvOverlap = true;
        }

        // proceed if there are 2+ out-of-frame INDELs or at least an INDEL and an overlapping SNV/MNV
        if(!hasInframeSnvMnvOverlap && frameshiftCount < 2)
            return;

        if(!isCodonMultiple(indelBaseTotal))
        {
            PV_LOGGER.trace("lps({}) varCount({}) out-of-frame with indel base count({})",
                    localPhaseSet, variants.size(), indelBaseTotal);
            return;
        }

        final String chromosome = variants.get(0).Chromosome;
        boolean posStrand = transImpacts.get(0).TransData.posStrand();

        String combinedRefCodons = "";
        String combinedAltCodons = "";
        int minCodonIndex = -1;
        int lastImpactedRefCodonEnd = 0;
        int lastRefCodonEnd = 0;

        for(int i = 0; i < transImpacts.size(); ++i)
        {
            VariantData variant = variants.get(i);

            if(ignoredVariants.contains(variant))
                continue;

            VariantTransImpact transImpact = transImpacts.get(i);

            if(!transImpact.hasCodingBases())
                continue;

            // a distinction is made between the recorded ref codon bases and those actually involved in / impact by the variant
            // for the purposes of checking for an overlap, the impacted bases are compared
            ImpactedRefCodingData transRefCodingData = getImpactedRefCodingData(refCodingImpacts, transImpact);
            int impactedRefCodonStart = transRefCodingData.PosStart;
            int impactedRefCodonEnd = transRefCodingData.PosEnd;
            int refCodonStart = transImpact.proteinContext().refCodingBaseStart();
            int refCodonEnd = transImpact.proteinContext().refCodingBaseEnd();

            VariantTransImpact nextTransImpact = i < transImpacts.size() - 1 ? transImpacts.get(i + 1) : null;

            int nextImpactedRefCodonStart = 0;

            if(nextTransImpact != null && nextTransImpact.hasCodingBases())
            {
                ImpactedRefCodingData nextRefCodingData = getImpactedRefCodingData(refCodingImpacts, nextTransImpact);
                nextImpactedRefCodonStart = nextRefCodingData.PosStart;
            }

            boolean overlapsOnStart = lastImpactedRefCodonEnd > 0 && impactedRefCodonStart <= lastImpactedRefCodonEnd;
            boolean overlapsOnEnd = nextImpactedRefCodonStart > 0 && impactedRefCodonEnd >= nextImpactedRefCodonStart;

            PV_LOGGER.trace("lps({}) var({}) codons({} -> {}) range({} - {}) overlaps(start={} end={})",
                    localPhaseSet, variant, transImpact.proteinContext().RefCodonBases, transImpact.proteinContext().AltCodonBases,
                    impactedRefCodonStart, impactedRefCodonEnd, overlapsOnStart, overlapsOnEnd);

            if(minCodonIndex < 0)
                minCodonIndex = transImpact.proteinContext().CodonIndex;
            else
                minCodonIndex = min(transImpact.proteinContext().CodonIndex, minCodonIndex);

            if(overlapsOnStart)
            {
                // first build the ref codons
                int refOverlap = lastRefCodonEnd - refCodonStart + 1;

                if(refOverlap <= transImpact.proteinContext().RefCodonBases.length())
                    combinedRefCodons += transImpact.proteinContext().RefCodonBases.substring(refOverlap);

                // add the alt codons, which do not extend any further position-wise than the ref codons and also end at a codon boundary
                // take the alt bases from the current since they will contain the change

                int lowerNonAltBases = variant.Position - refCodonStart;

                int currentAltBasesTrimmed;
                int prevAltBasesTrimmed;

                if(refOverlap >= lowerNonAltBases)
                {
                    currentAltBasesTrimmed = lowerNonAltBases;
                    prevAltBasesTrimmed = refOverlap - lowerNonAltBases;
                }
                else
                {
                    currentAltBasesTrimmed = refOverlap;
                    prevAltBasesTrimmed = 0;
                }

                if(combinedAltCodons.length() - prevAltBasesTrimmed <= 0)
                {
                    PV_LOGGER.warn("phasing variants LPS({}) var({}) combinedAltCodons({}) prevAltBasesTrimmed({})",
                            localPhaseSet, variant, combinedAltCodons, prevAltBasesTrimmed);
                    return;
                }

                String previousExtraAltBases = "";

                if(lastRefCodonEnd > refCodonEnd)
                {
                    previousExtraAltBases = combinedAltCodons.substring(lastRefCodonEnd - refCodonEnd + 1);
                }

                combinedAltCodons = combinedAltCodons.substring(0, combinedAltCodons.length() - prevAltBasesTrimmed);

                combinedAltCodons += transImpact.proteinContext().AltCodonBases.substring(currentAltBasesTrimmed);

                // restore the trimmed ref bases which the previous variant(s) had
                combinedAltCodons += previousExtraAltBases;
            }
            else
            {
                // fill in any missing gaps in the codons
                if(lastRefCodonEnd > 0 && refCodonStart > lastRefCodonEnd + 1)
                {
                    String gapCodonBases = refGenome.getBaseString(chromosome, lastRefCodonEnd + 1, refCodonStart - 1);
                    combinedRefCodons += gapCodonBases;
                    combinedAltCodons += gapCodonBases;
                }

                combinedRefCodons += transImpact.proteinContext().RefCodonBases;
                combinedAltCodons += transImpact.proteinContext().AltCodonBases;
            }

            lastRefCodonEnd = max(refCodonEnd, lastRefCodonEnd);
            lastImpactedRefCodonEnd = max(impactedRefCodonEnd, lastImpactedRefCodonEnd);
        }

        if(!isCodonMultiple(combinedRefCodons.length()) || !isCodonMultiple(combinedAltCodons.length()))
            return;

        ProteinContext combinedPc = new ProteinContext();
        combinedPc.CodonIndex = minCodonIndex;
        combinedPc.RefCodonBases = combinedRefCodons;
        combinedPc.AltCodonBases = combinedAltCodons;

        // check synonymous vs missense from combined bases
        VariantEffect combinedEffect;

        if(posStrand)
        {
            combinedPc.RefAminoAcids = Codons.aminoAcidFromBases(combinedRefCodons);
            combinedPc.AltAminoAcids = Codons.aminoAcidFromBases(combinedAltCodons);
        }
        else
        {
            combinedPc.RefAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(combinedRefCodons));
            combinedPc.AltAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(combinedAltCodons));
        }

        trimAminoAcids(combinedPc, true, true, false);

        if(indelBaseTotal == 0)
        {
            combinedEffect = combinedPc.RefAminoAcids.equals(combinedPc.AltAminoAcids) ? PHASED_SYNONYMOUS : PHASED_MISSENSE;

            PV_LOGGER.trace("lps({}) varCount({}) combinedEffect({}) from aminoAcids({} -> {})",
                    localPhaseSet, variants.size(), combinedEffect, combinedPc.RefAminoAcids, combinedPc.AltAminoAcids);
        }
        else
        {
            combinedEffect = indelBaseTotal > 0 ? PHASED_INFRAME_INSERTION : PHASED_INFRAME_DELETION;

            PV_LOGGER.trace("lps({}) varCount({}) indelBaseTotal({})",
                    localPhaseSet, variants.size(), indelBaseTotal);
        }

        List<VariantEffect> combinedEffects = Lists.newArrayList(combinedEffect);

        VariantEffect ssEffect = checkStopStartCodons(minCodonIndex, combinedPc.RefAminoAcids, combinedPc.AltAminoAcids);

        if(ssEffect != null)
        {
            if(ssEffect == STOP_LOST || ssEffect == START_LOST)
            {
                combinedEffects.remove(MISSENSE); // superceded if present so remove
                combinedEffects.remove(FRAMESHIFT); // superceded if present so remove
            }

            combinedEffects.add(ssEffect);
        }

        combinedPc.IsPhased = true;
        combinedPc.Hgvs = HgvsProtein.generate(combinedPc, combinedEffects);

        // now convert missense / synonymous to phased inframe so it's clearer what has happened
        if(combinedEffects.contains(SYNONYMOUS) || combinedEffects.contains(MISSENSE))
        {
            combinedEffects.remove(MISSENSE);
            combinedEffects.remove(SYNONYMOUS);
            combinedEffects.add(indelBaseTotal > 0 ? PHASED_INFRAME_INSERTION : PHASED_INFRAME_DELETION);
        }

        for(int i = 0; i < transImpacts.size(); ++i)
        {
            VariantData variant = variants.get(i);

            // ignore SNVs/MNVs if they don't overlap the same codons as any of the INDELs unless they are between indels
            if(ignoredVariants.contains(variant))
                continue;

            VariantTransImpact transImpact = transImpacts.get(i);

            if(!transImpact.hasCodingBases())
                continue;

            transImpact.markPhasedFrameshift();
            transImpact.setProteinContext(combinedPc);

            transImpact.codingContext().IsFrameShift = false;

            transImpact.effects().clear();
            combinedEffects.forEach(x -> transImpact.addEffect(x));
        }
    }

    private class ImpactedRefCodingData
    {
        public final int PosStart;
        public final int PosEnd;

        public ImpactedRefCodingData(final int posStart, final int posEnd)
        {
            PosStart = posStart;
            PosEnd = posEnd;
        }

        public String toString() { return format("%d - %d", PosStart, PosEnd); }
    }

    private ImpactedRefCodingData getImpactedRefCodingData(
            final Map<VariantTransImpact,ImpactedRefCodingData> impactMap, final VariantTransImpact transImpact)
    {
        ImpactedRefCodingData impactData = impactMap.get(transImpact);

        if(impactData == null)
        {
            impactData = impactedRefCodingBasePosition(transImpact);
            impactMap.put(transImpact, impactData);
        }

        return impactData;
    }

    private ImpactedRefCodingData impactedRefCodingBasePosition(final VariantTransImpact transImpact)
    {
        boolean inframeIndel = transImpact.hasEffect(INFRAME_DELETION) || transImpact.hasEffect(INFRAME_INSERTION);
        final ProteinContext pc = transImpact.proteinContext();

        int refCodingBaseStart = pc.refCodingBasePosition(SE_START);
        int refCodingBaseEnd = pc.refCodingBasePosition(SE_END);

        if(!inframeIndel)
        {
            return new ImpactedRefCodingData(refCodingBaseStart, refCodingBaseEnd);
        }

        // special case for inframe INDELs which have included a full ref codon at the start and/or at end to aid with AA evaluation
        int netCodonImpact = pc.NetCodonIndexRange[SE_END] - pc.NetCodonIndexRange[SE_START] + 1;
        int refCodonLengthDiff = abs(pc.RefCodonBases.length() - pc.AltCodonBases.length());

        if(netCodonImpact > 0 && refCodonLengthDiff == netCodonImpact * CODON_LENGTH)
        {
            // strip the unaffected ref codon from the start and/or end
            if(netCodonImpact > 1)
            {
                refCodingBaseStart += CODON_LENGTH;
                refCodingBaseEnd -= CODON_LENGTH;
            }
            else
            {
                if(refCodingBaseStart == transImpact.codingContext().CodingPositionRange[SE_START])
                    refCodingBaseEnd -= CODON_LENGTH;
                else
                    refCodingBaseStart += CODON_LENGTH;
            }
        }

        // no extension was required so use the standard positions
        return new ImpactedRefCodingData(refCodingBaseStart, refCodingBaseEnd);
    }
}
