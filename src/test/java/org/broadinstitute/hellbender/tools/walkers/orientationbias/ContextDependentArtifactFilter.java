package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Map;
import java.util.TreeMap;

/**
 * Created by tsato on 7/26/17.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = VariantProgramGroup.class
)
public class ContextDependentArtifactFilter extends VariantWalker {
    // Consider using an array instead; then we need a contextToIndex() method
    final Map<String, Integer> contextCountMap = new TreeMap();
    final int numVariants = 300; // any way to learn this a priori?
    final int[] altCount = new int[numVariants];
    // Isn't there a context data structure - maybe I saw this in Lee's code
    final String[] context = new String[numVariants];


    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        String tumorSampleName = "synthetic.challenge.set1.tumor"; // TODO: get it from vc or vcf
        String referenceBases = variant.getAttributeAsString(GATKVCFConstants.REFERENCE_CONTEXT_KEY, "");
        if (referenceBases.equals("")){
            // log an error and move on
        }

        final int ALT_INDEX = 1;
        int altReadCount = variant.getGenotype(tumorSampleName).getAD()[ALT_INDEX]; // alternate allele

        // what are the sufficient statistics?

    }

    @Override
    public Object onTraversalSuccess(){

        return null;
    }
}
