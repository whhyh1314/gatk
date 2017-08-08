package org.broadinstitute.hellbender.tools.walkers.orientationbias;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.IndexRange;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Created by tsato on 7/26/17.
 */
public class ContextDependentArtifactFilterEngine {
    // z \in { F1R2, F2R1, Balanced Hom Ref, Balanced Het, Balanced Hom Var }. Thus |z| = 5.
    static final int NUM_STATUSES = 5;

    // 4 possible bases in 3 positions, so 4^3 = 64
    static final int NUM_CONTEXTS = 64;

    // A, C, G, or T, since we only look at SNP sites
    static final int NUM_POSSIBLE_ALLELES = 4;

    // hyperparameters of the model
    // pi is the weight vector (prior?) for the categorical variable z.
    // we have a NUM_STATUSES-dimensional vector for each of
    // TODO: consider using a matrix here
    double[][][] pi = new double[NUM_CONTEXTS][NUM_POSSIBLE_ALLELES][NUM_STATUSES];

    // probability of heads for the binomial m, which represents the number of alt reads at site n
    double[] f = new double[NUM_CONTEXTS];

    // probability of heads for the binomial x, which represents the number of F1R2 alt reads at site n
    double[] theta = new double[NUM_CONTEXTS];

    // number of 3-mers over which we have reads in the input bam
    int numLoci;

    // the posterior probability of latent variable z with the current estimates of hyperparameters pi, f and theta
    RealMatrix responsibilities;

    int numIterations = 0;

    // When the increase in likelihood falls below this value we deem the algorithm converged
    final double EPSILON = 1e-4;

    public ContextDependentArtifactFilterEngine(){
        initializeParameters();
    }

    public void runEMAlgorithm(final short[] observedAltReadCounts, final short[] observedAltF1R2Counts,
                               final int[] observedDepths, final String[] observedRefContexts){
        // remember we run EM separately for each of 4^3 = 64 ref contexts
        List<String> refContexts = new ArrayList<>();

        // TODO: systematically create all reference contexts
        refContexts.add(0, "ATA");
        refContexts.add(1, "ACG");
        for (String refContext : refContexts){
            // only use the subset of data points with this particular ref context
            final int[] indices = IntStream.range(0, numLoci).filter(i -> observedRefContexts[i].equals(refContexts)).toArray();
            final int[] depths = Arrays.stream(indices).map(i -> observedDepths[i]).toArray();
            final int[] altReadCounts = Arrays.stream(indices).map(i -> observedAltReadCounts[i]).toArray();
            final int[] altF1R2Counts = Arrays.stream(indices).map(i -> observedAltF1R2Counts[i]).toArray();

            while (! checkLikelihoodHasConverged(numIterations)) {
                takeMstep();

                takeEstep();
                numIterations++;
            }

        }



    }

    private void initializeParameters(){
        responsibilities = new Array2DRowRealMatrix(numLoci, NUM_STATUSES);
        // initialize responsibilities to be flat
        final double initialProbability = 1.0 / NUM_STATUSES;
        for (int i = 0; i < numLoci; i++ ) {
            for (int j = 0; j < NUM_STATUSES; j++) {
                responsibilities.setEntry(i, j, initialProbability);
            }
        }


    }

    private void takeEstep(){

    }

    private void takeMstep(){
        // TODO: should it take responsibilities as an argument?


    }

    private boolean checkLikelihoodHasConverged(final double oldLikelihood, final double newLikelihood){
        return Math.abs(newLikelihood - oldLikelihood) < EPSILON;
    }

    private boolean checkLikelihoodHasConverged(final int numIterations){
        return numIterations > 10;
    }

}
