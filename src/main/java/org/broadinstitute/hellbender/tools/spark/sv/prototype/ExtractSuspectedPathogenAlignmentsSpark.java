package org.broadinstitute.hellbender.tools.spark.sv.prototype;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

/**
 * Tool to select alignments of long reads possibly pointing to pathogen integration site.
 */
@CommandLineProgramProperties(summary="Tool to select alignments of long reads possibly pointing to pathogen integration site.",
        oneLineSummary="Tool to select alignments of long reads possibly pointing to pathogen integration site.",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public class ExtractSuspectedPathogenAlignmentsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(ExtractSuspectedPathogenAlignmentsSpark.class);

    @Argument(doc = "sam file for aligned contigs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputAssemblyAlignments;

    @Argument(doc = "length of clip a uncovered read must have for it to be included in output", shortName = "uci",
            fullName = "uncoveredClipLength")
    private int uncoveredClipLength;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }


    @Override
    protected void runTool(final JavaSparkContext ctx) {

        select(ctx, getUnfilteredReads(), uncoveredClipLength, getHeaderForReads(), outputAssemblyAlignments, localLogger);
    }

    JavaRDD<AlignedContig> select(final JavaSparkContext ctx, final JavaRDD<GATKRead> rawAlignments,
                                  final int coverageThresholdInclusive, final SAMFileHeader header, final String selectedAlignmentFileName,
                                  final Logger toolLogger) {
        final JavaRDD<GATKRead> unmappedReads = rawAlignments.filter(GATKRead::isUnmapped);

        final JavaRDD<AlignedContig> selectedMappedContigs =
                new DiscoverVariantsFromContigAlignmentsSAMSpark
                        .SAMFormattedContigAlignmentParser(rawAlignments.filter(r -> !r.isUnmapped()), header, true, localLogger)
                        .getAlignedContigs()
                        .filter(ctg -> keepContigForPathSeqUse(ctg, coverageThresholdInclusive))
                        .cache();

        final HashSet<String> readNames = new HashSet<>(selectedMappedContigs.map(tig -> tig.contigName).distinct().collect());
        final JavaRDD<GATKRead> readsToWrite = rawAlignments.filter(r -> readNames.contains(r.getName())).union(unmappedReads);
        if (readsToWrite.isEmpty()) toolLogger.warn("Found no alignments indicating possible pathogen integration sites.");
        writeSAMFile(readsToWrite.map(r -> r.convertToSAMRecord(header)).toLocalIterator(), header, selectedAlignmentFileName, false);

        return selectedMappedContigs;
    }

    /**
     * Currently keep contigs that has
     *   1) alignments covering less than or equal to only half of the contig sequence,
     *   2) a certain length of its sequence uncovered by alignments
     */
    private static boolean keepContigForPathSeqUse(final AlignedContig contig,
                                                   final int coverageThresholdInclusive) {

        final int alignmentCoverage = alignmentsCoverage(contig);

        return contig.contigSequence.length >= 2 * alignmentCoverage ||
                contig.contigSequence.length - alignmentCoverage >= coverageThresholdInclusive;
    }

    /**
     * Computes how long of the provided contig is covered with alignments.
     */
    @VisibleForTesting
    static int alignmentsCoverage(final AlignedContig contig) {
        if (contig.alignmentIntervals.isEmpty()) return 0;

        final List<AlignmentInterval> alignmentIntervals = contig.alignmentIntervals;
        final List<SVInterval> maximallyExtendedCovers = new ArrayList<>(alignmentIntervals.size());

        final Iterator<AlignmentInterval> it = alignmentIntervals.iterator();
        AlignmentInterval current = it.next(); // at least two exist
        SVInterval currentSVI = new SVInterval(1, current.startInAssembledContig, current.endInAssembledContig+1); // +1 for [a,b) semi-closed
        while(it.hasNext()) {
            AlignmentInterval next = it.next();
            final SVInterval nextSVI = new SVInterval(1, next.startInAssembledContig, next.endInAssembledContig+1);

            if (nextSVI.overlaps(currentSVI)) {
                currentSVI = new SVInterval(1, Math.min(currentSVI.getStart(), nextSVI.getStart()),
                        Math.max(currentSVI.getEnd(), nextSVI.getEnd())     );
            } else {
                maximallyExtendedCovers.add(currentSVI);
                currentSVI = nextSVI;
            }
        }
        maximallyExtendedCovers.add(currentSVI);

        return maximallyExtendedCovers.stream().mapToInt(SVInterval::getLength).sum();
    }

    public static void writeSAMFile(final Iterator<SAMRecord> alignments, final SAMFileHeader header, final String outputName,
                                    final boolean preOrdered) {
        try ( final SAMFileWriter writer = createSAMFileWriter(outputName, header, preOrdered) ) {
            alignments.forEachRemaining(writer::addAlignment);
        } catch ( final UncheckedIOException ie) {
            throw new GATKException("Can't write SAM file to the specified location: " + outputName, ie);
        }
    }

    private static SAMFileWriter createSAMFileWriter(final String samFile, final SAMFileHeader header, final boolean preOrdered) {
        final SAMFileWriterFactory factory = new SAMFileWriterFactory();
        final int lastDotIndex = samFile.lastIndexOf('.');
        if (lastDotIndex >= 0) {
            final String extension = samFile.substring(lastDotIndex).toLowerCase();
            if (extension.equals(BamFileIoUtils.BAM_FILE_EXTENSION)) {
                return factory.makeBAMWriter(header, preOrdered, BucketUtils.createFile(samFile));
            } else if (extension.equals(".sam")) {
                return factory.makeSAMWriter(header, preOrdered, BucketUtils.createFile(samFile));
            } else {
                throw new GATKException("unsupported read alignment file name extension (." + extension + ") in requested name: " + samFile);
            }
        } else {
            throw new GATKException("cannot determine the alignment file format from its name: " + samFile);
        }
    }
}