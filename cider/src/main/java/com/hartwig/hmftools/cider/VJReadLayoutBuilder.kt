package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.LayoutTree
import com.hartwig.hmftools.cider.layout.ReadLayout
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.ArrayList

// helper class to convert from the outer VJ classes to the layout classes
// create an interface to make it easier to test
abstract class IVJReadLayoutAdaptor
{
    abstract fun toReadCandidate(read: ReadLayout.Read) : VJReadCandidate

    abstract fun toLayoutReadSlice(read: ReadLayout.Read) : ReadSlice

    abstract fun getAnchorMatchMethod(layout: ReadLayout) : VJReadCandidate.MatchMethod
    abstract fun getTemplateAnchorSequence(layout: ReadLayout) : String

    // anchor range that can be outside the layout, by extrapolating it to what
    // position the anchor should be
    abstract fun getExtrapolatedAnchorRange(vj: VJ, layout: ReadLayout) : IntRange

    // returns the part of the anchor range that is inside the layout, null if it is not
    // inside
    fun getAnchorRange(vj: VJ, layout: ReadLayout) : IntRange?
    {
        val anchorRange = getExtrapolatedAnchorRange(vj, layout)

        if (anchorRange.first >= layout.length || anchorRange.last < 0)
            return null

        // protect against 0 and end
        return Math.max(0, anchorRange.first) until Math.min(anchorRange.last + 1, layout.length)
    }
}

// This class is the mapper between the candidate reads and the layout object
// There are several hacks that is used here to keep track of
//
// 1. for layout that comes from V read candidates, the aligned position is the last base of the V anchor i.e.
//                                  * <-- this T is the aligned position of the layout
//    AGATCTGAG-GACACGGCCGTGTATTACTGT-GCGAGAGACACAGTGTGAAAACCCACATCCTGAGAGTGTCAGAAACCCTGAGGGA
//              |___________________|
//                V anchor
// 2. for layout that comes from J read candidates, the aligned position is the first base of the J anchor, i.e.
//         this C is the aligned position of the layout --> *
//    AGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGACACAGTGTGAAAACC-CACATCCTGAGAGTGTCAGAA-ACCCTGAGGGA
//                                                          |___________________|
//                                                                 J anchor
// Most functions here rely on this.
class VJReadLayoutBuilder(private val trimBases: Int, private val minBaseQuality: Int) : IVJReadLayoutAdaptor()
{
    private class VjLayoutRead private constructor(
        val layoutReadSlice: ReadSlice,
        val readCandidate: VJReadCandidate,
        readKey: ReadKey,
        sequence: String,
        baseQualities: ByteArray,
        alignedPosition: Int)
        : ReadLayout.Read(readKey, sequence, baseQualities, alignedPosition)
    {
        constructor(layoutReadSlice: ReadSlice, readCandidate: VJReadCandidate, alignedPosition: Int)
                : this(layoutReadSlice,
            readCandidate,
            ReadKey(layoutReadSlice.readName, layoutReadSlice.firstOfPairFlag),
            layoutReadSlice.readString,
            layoutReadSlice.baseQualities,
            alignedPosition)
        {
        }

        override fun copy(alignedPosition: Int): ReadLayout.Read
        {
            return VjLayoutRead(layoutReadSlice, readCandidate, readKey, sequence, baseQualities, alignedPosition)
        }
    }

    fun readCandidateToLayoutRead(readCandidate: VJReadCandidate) : ReadLayout.Read?
    {
        val slice = determineReadSlice(readCandidate.read, readCandidate.useReverseComplement)

        if (slice == null)
            return null

        // now determine the aligned position
        val alignedPosition: Int

        if (readCandidate.vjGeneType.vj == VJ.V)
        {
            // aligned position we must take into account that we remove all bases before vAnchor start
            alignedPosition = readCandidate.anchorOffsetEnd - 1 - slice.sliceStart
        }
        else
        {
            // aligned position we must take into account that we remove all bases before vAnchor start
            alignedPosition = readCandidate.anchorOffsetStart - slice.sliceStart
        }

        return VjLayoutRead(slice, readCandidate, alignedPosition)
    }

    override fun toReadCandidate(read: ReadLayout.Read) : VJReadCandidate
    {
        return (read as VjLayoutRead).readCandidate
    }

    override fun toLayoutReadSlice(read: ReadLayout.Read) : ReadSlice
    {
        return (read as VjLayoutRead).layoutReadSlice
    }

    fun getReadCandidates(layout: ReadLayout) : List<VJReadCandidate>
    {
        return layout.reads.map({ read: ReadLayout.Read -> toReadCandidate(read)})
    }

    // get the anchor boundary position for this layout read
    fun getAnchorBoundaryPosition(read: ReadLayout.Read) : Int
    {
        return read.alignedPosition
    }

    // apply trim bases and polyG trimming
    private fun determineReadSlice(read: SAMRecord, useReverseComplement: Boolean) : ReadSlice?
    {
        // work out the slice start and end
        var sliceStart: Int = trimBases
        var sliceEnd: Int = read.readLength - trimBases

        // now we also want to try poly G tail trimming
        // we want to work out there the tail is.
        // the tail is on the right side and poly G if !read.readNegativeStrandFlag
        // the tail is on the left side and poly C otherwise
        if (!read.readNegativeStrandFlag)
        {
            // ends with poly G, but take trim bases into account
            val numGs = CiderUtils.numTrailingPolyG(read.readString, sliceEnd)
            if (numGs >= CiderConstants.MIN_POLY_G_TRIM_COUNT)
            {
                sLogger.debug("read({}) strand(+) poly G tail of length({}) found({})",
                    read, numGs, read.readString)
                sliceEnd -= numGs + CiderConstants.POLY_G_TRIM_EXTRA_BASE_COUNT
            }
        }
        else
        {
            val numCs = CiderUtils.numLeadingPolyC(read.readString, sliceStart)
            if (numCs >= CiderConstants.MIN_POLY_G_TRIM_COUNT)
            {
                sLogger.debug("read({}) strand(-) poly G tail of length({}) found({})",
                    read, numCs, read.readString)
                sliceStart += numCs + CiderConstants.POLY_G_TRIM_EXTRA_BASE_COUNT
            }
        }

        // the above logic is before reverse complement, the following logic is after
        // so we swap the start / end here
        if (useReverseComplement)
        {
            val sliceStartTmp = sliceStart
            sliceStart = read.readLength - sliceEnd
            sliceEnd = read.readLength - sliceStartTmp
        }

        if ((sliceEnd - sliceStart) < 5)
        {
            // if too little left don't bother
            return null
        }

        return ReadSlice(read, useReverseComplement, sliceStart, sliceEnd)
    }

    override fun getAnchorMatchMethod(layout: ReadLayout): VJReadCandidate.MatchMethod
    {
        val readCandidates = getReadCandidates(layout)

        // we need to get a few values from read candidates
        if (readCandidates.isEmpty())
            throw IllegalArgumentException("read candidate list is empty")

        // just return first one for now, not the best but should be fine
        return readCandidates.first().matchMethod
    }

    override fun getTemplateAnchorSequence(layout: ReadLayout) : String
    {
        val readCandidates = getReadCandidates(layout)

        // we need to get a few values from read candidates
        if (readCandidates.isEmpty())
            throw IllegalArgumentException("read candidate list is empty")

        return readCandidates.maxByOrNull({ r -> r.similarityScore })?.templateAnchorSequence ?: ""
    }

    // NOTE: the anchor range could be outside the layout
    override fun getExtrapolatedAnchorRange(vj: VJ, layout: ReadLayout) : IntRange
    {
        val layoutReads = layout.reads.map { o: ReadLayout.Read -> toReadCandidate(o) }
            .toList()

        // we more or less get the top one
        val anchorLength = layoutReads.maxOfOrNull { o: VJReadCandidate -> o.anchorOffsetEnd - o.anchorOffsetStart } ?: 0

        val anchorRange =
        // for V read we align to last base of anchor, for J read we align to first base of the anchor
        when (vj)
        {
            VJ.V -> layout.alignedPosition - anchorLength + 1..layout.alignedPosition
            VJ.J -> layout.alignedPosition until layout.alignedPosition + anchorLength
        }

        return anchorRange
    }

    // Build layouts using layout tree
    // TODO: merge layouts where differences are outside of the anchor
    fun buildLayouts(geneType: VJGeneType, readCandidates: List<VJReadCandidate>,
                     minMatchedBases: Int)
    : List<ReadLayout>
    {
        sLogger.info("building {} layouts from {} reads", geneType, readCandidates.size)

        val layoutReads = ArrayList<ReadLayout.Read>()

        for (r in readCandidates)
        {
            val layoutRead = readCandidateToLayoutRead(r)
            if (layoutRead != null)
                layoutReads.add(layoutRead)
        }

        // always build from left to right
        layoutReads.sortWith(
            Collections.reverseOrder(
                Comparator.comparingInt({ r: ReadLayout.Read -> r.alignedPosition })
                    .thenComparingDouble({ r: ReadLayout.Read -> r.baseQualities.average() }) // handle the highest quality ones first
                    .thenComparingInt({ r: ReadLayout.Read -> r.alignedPosition + r.sequence.length })
                    .thenComparing({ r: ReadLayout.Read -> r.readKey.readName }) // lastly we use read Id just in case
            ))

        sLogger.info("building layouts from {} reads", layoutReads.size)

        val layoutTree = LayoutTree(minBaseQuality.toByte(), minMatchedBases)

        // go through the read data list, and add one by one to the list of clusters
        // if there are multiple clusters that matches, we choose the highest one
        for (read in layoutReads)
        {
            layoutTree.tryAddRead(LayoutTree.Read(read, read.sequence, read.baseQualities, read.alignedPosition))
        }

        val readLayouts: List<ReadLayout> = layoutTree.buildReadLayouts({ read: LayoutTree.Read -> read.source as VjLayoutRead })

        sLogger.info("built {} layouts from {} reads", readLayouts.size, layoutReads.size)

        return readLayouts
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(VJReadLayoutBuilder::class.java)
    }
}