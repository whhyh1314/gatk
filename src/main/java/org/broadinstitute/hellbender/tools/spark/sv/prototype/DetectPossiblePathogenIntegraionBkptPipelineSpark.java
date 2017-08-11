package org.broadinstitute.hellbender.tools.spark.sv.prototype;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.FindBreakpointEvidenceSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Tool to run the sv pipeline up for possible pathogen integration site detection and assembled contigs with low alignment coverage.
 */
@CommandLineProgramProperties(summary="Tool to run the sv pipeline up for possible pathogen integration site detection and assembled contigs with low alignment coverage.",
        oneLineSummary="Tool to run the sv pipeline up for possible pathogen integration site detection and assembled contigs with low alignment coverage.",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public class DetectPossiblePathogenIntegraionBkptPipelineSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(DetectPossiblePathogenIntegraionBkptPipelineSpark.class);

    @Argument(doc = "sam file for aligned contigs", shortName = "contigSAMFile",
            fullName = "contigSAMFile")
    private String outputAssemblyAlignments;

    @Argument(doc = "filename for output vcf", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutputFileName;


    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection evidenceAndAssemblyArgs
            = new StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection();

    @ArgumentCollection
    private StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs
            = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

    @Argument(doc = "length of clip a uncovered read must have for it to be included in output", shortName = "uci",
            fullName = "uncoveredClipLength")
    private int uncoveredClipLength;

    @Argument(doc = "whether to write alignments of all successful assemblies", shortName = "allAln",
            fullName = "allAln")
    private boolean writeAllAlignments = false;

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

        final SAMFileHeader header = getHeaderForReads();
        final PipelineOptions pipelineOptions = getAuthenticatedGCSOptions();

        final List<String> refNames = AlignedAssemblyOrExcuse.getRefNames(header);

        // gather evidence, run assembly, and align, but NOT outputting the alignments yet unless specifically asked to
        final JavaRDD<GATKRead> rawAlignments =
                ctx.parallelize(
                FindBreakpointEvidenceSpark
                        .gatherEvidenceAndWriteContigSamFile(ctx, evidenceAndAssemblyArgs, header, getUnfilteredReads(),
                                writeAllAlignments ? outputAssemblyAlignments.replace(".sam", ".all.sam").replace(".bam", ".all.bam") : null,
                                localLogger).stream()
                        .filter(AlignedAssemblyOrExcuse::isNotFailure)
                        .flatMap(aaoe -> aaoe.toSAMStreamForAlignmentsOfThisAssembly(header, refNames))
                        .map(SAMRecordToGATKReadAdapter::new)
                        .collect(Collectors.toList()));
        if (rawAlignments.isEmpty()) return;

        // parse the contig alignments and select alignments
        final int coverageThresholdInclusive = uncoveredClipLength;
        final JavaRDD<AlignedContig> selectedAlignments =
                new ExtractSuspectedPathogenAlignmentsSpark().select(ctx, rawAlignments, coverageThresholdInclusive,
                        header, outputAssemblyAlignments, localLogger);
        if(selectedAlignments.isEmpty()) {
            localLogger.warn("Could not find any sites suggesting pathogen integration.");
            return;
        }

        // discover variants and write to vcf
        DiscoverVariantsFromContigAlignmentsSAMSpark
                .discoverVariantsAndWriteVCF(selectedAlignments, discoverStageArgs.fastaReference,
                        ctx.broadcast(getReference()), pipelineOptions, vcfOutputFileName, localLogger);
    }

}
