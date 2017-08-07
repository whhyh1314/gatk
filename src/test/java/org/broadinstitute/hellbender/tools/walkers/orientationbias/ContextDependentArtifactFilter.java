package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.List;
import java.util.Optional;

/**
 * Created by tsato on 7/26/17.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = VariantProgramGroup.class
)
public class ContextDependentArtifactFilter extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "", common = false, optional = true)
    private File OUTPUT_FILE = null;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "", optional = true)
    private FeatureInput<VariantContext> mutect2VcfFile;

    @Argument(fullName = "", shortName = "", doc = "", optional = true)
    static int DEFAULT_INITIAL_LIST_SIZE = 37_000_000/64; // by default we assume that that all 64 reference 3-mers are equally likely

    @Override
    public boolean requiresReference(){
        return true;
    }




    // Assume that we have a whole exome bam (~37M bases) plus give it a 1M base cushion
    final int numVariants = 37_000_000 + 1_000_000;

    /**
     * Arrays of observed variables and parameters. We have > 1 million entries so we should pick the data type of
     * the arrays carefully. The max of short primitive type is Short.MAX_VALUE = 32,767,
     * which should be plenty for observedAltReadCounts and altF1R2count. We will use int[] for depth to be
     * on the safe side
     **/

    // TODO: note that most of these counts are zero - think of a more sparse encoding here
    // m in the graphical model
    final short[] observedAltReadCounts = new short[numVariants];

    // x in the graphical model
    final short[] observedAltF1R2Counts = new short[numVariants];

    // R in the graphical model
    final int[] observedDepths = new int[numVariants];

    // 3-mer reference context - TODO: use Kmer? Is it more efficient to create an array of object, so that the
    // size of an object is predictable, than to make string, for which we store pointers which will point to
    // places all over the memory?
    final String[] observedRefContexts = new String[numVariants];

    // n in [0, N), where N is the number of mutect2VcfFile in the vcf
    int n = 0;

    private ContextDependentArtifactFilterEngine engine;

    @Override
    public void onTraversalStart(){
        FeatureDataSource<VariantContext> vcfFile = new FeatureDataSource(new File(mutect2VcfFile.getFeaturePath()));
    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext){
        // Does it only traverse over vcfs this way?
        final ReadPileup pileup = alignmentContext.getBasePileup();
        final int start = alignmentContext.getStart();
        List<VariantContext> variantContexts = featureContext.getValues(mutect2VcfFile);
        SimpleInterval window = referenceContext.getWindow();

        // referenceContext always comes withe window of single base, so
        // manually expand the window and get the 3-mer for now.
        // TODO: implement getBasesInInterval() in referenceContext. Maybe simplify to getKmer(int k)?
        referenceContext.setWindow(1, 1);
        final String referenceKmer = new String(referenceContext.getBases());

        final int[] baseCounts = pileup.getBaseCounts();

        // R in the docs
        final int depth = (int) MathUtils.sum(baseCounts);
        observedDepths[n] = depth;

        observedRefContexts[n] = referenceKmer;

        if (variantContexts.isEmpty()){
            // no variants at this locus
            observedAltReadCounts[n] = 0;
            observedAltF1R2Counts[n] = 0;
            n++;
            return;
        }

        // We have a variant at this locus, and we assume that we always get one variant context.
        // Also assume that the site is single-allelic
        final Optional<Allele> altAllele = variantContexts.isEmpty() ? Optional.empty() :
                Optional.of(variantContexts.get(0).getAlternateAllele(0));

        if (altAllele.get().getBases().length > 1){
            // How should we handle INDEL sites?
            // TODO: how should we avoid other messy sites? Refer to David's CalculateContamination for how he skipped certain sites
            return;
        }

        final byte altBase = altAllele.get().getBases()[0];

        // m in the docs. We can safely assume here that the variant is a SNP
        observedAltReadCounts[n] = (short) baseCounts[BaseUtils.simpleBaseToBaseIndex(altBase)];
        observedAltF1R2Counts[n] = (short) pileup.getReads().stream().filter(r -> r.getBase(0) == altBase)
                .filter(r -> ReadUtils.isF2R1(r))
                .count();

        n++;
        return;
    }

    @Override
    public Object onTraversalSuccess() {
        final int numLocus = n;
        engine = new ContextDependentArtifactFilterEngine(numLocus);
        engine.runEMAlgorithm(observedAltReadCounts, observedAltF1R2Counts, observedDepths, observedRefContexts);
        return null;
    }



}
