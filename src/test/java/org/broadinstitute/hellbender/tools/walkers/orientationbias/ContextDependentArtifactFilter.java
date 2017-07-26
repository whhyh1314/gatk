package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

/**
 * Created by tsato on 7/26/17.
 */
public class ContextDependentArtifactFilter extends VariantWalker {
    

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

    }
}
