package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Provide simulated AlleleFractionData objects and stores the "truth" data from which they were created
 *
 * @author David Benjamin
 */
public final class AlleleFractionSimulatedData {
    public static final ReadCountCollection TRIVIAL_TARGETS = new ReadCountCollection(
            Collections.singletonList(new Target("target", new SimpleInterval("chr99", 999999999, 999999999))),
            Collections.singletonList("SAMPLE"),
            new Array2DRowRealMatrix(new double[][] {{1}}));

    private static final int MIN_HETS_PER_SEGMENT = 3;
    private static final int RANDOM_SEED = 13;
    private static final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private static PoissonDistribution makePoisson(final RandomGenerator rng, final double mean) {
        return new PoissonDistribution(rng, mean, PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
    }

    private final AlleleFractionState trueState;
    private final PhaseIndicators truePhases = new PhaseIndicators(new ArrayList<>());
    private final SegmentedGenome segmentedGenome;
    private final int numSegments;

    public AlleleFractionSimulatedData(final double averageHetsPerSegment, final int numSegments,
            final double averageDepth, final double biasMean, final double biasVariance, final double outlierProbability) {
        rng.setSeed(RANDOM_SEED);
        this.numSegments = numSegments;
        final AlleleFractionState.MinorFractions minorFractions = new AlleleFractionState.MinorFractions(numSegments);
        final List<AllelicCount> alleleCounts = new ArrayList<>();
        final List<SimpleInterval> segments = new ArrayList<>();

        final PoissonDistribution segmentLengthGenerator = makePoisson(rng, averageHetsPerSegment);
        final PoissonDistribution readDepthGenerator = makePoisson(rng, averageDepth);
        final UniformRealDistribution minorFractionGenerator = new UniformRealDistribution(rng, 0.0, 0.5);

        //translate to ApacheCommons' parametrization of the gamma distribution
        final double gammaShape = biasMean * biasMean / biasVariance;
        final double gammaScale = biasVariance / biasMean;
        final GammaDistribution biasGenerator = new GammaDistribution(rng, gammaShape, gammaScale);

        //put each segment on its own chromosome and sort by lexicographical order
        final List<String> chromosomes = IntStream.range(0, numSegments).mapToObj(Integer::toString).collect(Collectors.toList());
        Collections.sort(chromosomes);

        for (final String chromosome : chromosomes) {
            // calculate the range of het indices for this segment
            final int numHetsInSegment = Math.max(MIN_HETS_PER_SEGMENT, segmentLengthGenerator.sample());

            final double minorFraction = minorFractionGenerator.sample();
            minorFractions.add(minorFraction);

            //we will put all the hets in this segment/chromosome at loci 1, 2, 3 etc
            segments.add(new SimpleInterval(chromosome, 1, numHetsInSegment + 1));
            for (int het = 1; het < numHetsInSegment + 1; het++) {
                final double bias = biasGenerator.sample();

                //flip a coin to decide alt minor (alt fraction = minor fraction) or ref minor (alt fraction = 1 - minor fraction)
                final boolean isAltMinor = rng.nextDouble() < 0.5;
                final double altFraction =  isAltMinor ? minorFraction : 1 - minorFraction;

                //the probability of an alt read is the alt fraction modified by the bias or, in the case of an outlier, random
                final double pAlt;
                if (rng.nextDouble() < outlierProbability) {
                    truePhases.add(AlleleFractionIndicator.OUTLIER);
                    pAlt = rng.nextDouble();
                } else {
                    truePhases.add(isAltMinor ? AlleleFractionIndicator.ALT_MINOR : AlleleFractionIndicator.REF_MINOR);
                    pAlt = altFraction / (altFraction + (1 - altFraction) * bias);
                }

                final int numReads = readDepthGenerator.sample();
                final int numAltReads = new BinomialDistribution(rng, numReads, pAlt).sample();
                final int numRefReads = numReads - numAltReads;
                alleleCounts.add(new AllelicCount(new SimpleInterval(chromosome, het, het), numRefReads, numAltReads));
            }
        }

        final Genome genome = new Genome(TRIVIAL_TARGETS, alleleCounts);
        segmentedGenome = new SegmentedGenome(segments, genome);
        trueState = new AlleleFractionState(biasMean, biasVariance, outlierProbability, minorFractions);
    };


    public AlleleFractionState getTrueState() { return trueState; }

    /**
     * Returns the ArrayList of phase indicators held internally, which should not be modified by the caller.
     */
    public PhaseIndicators getTruePhases() {
        return truePhases;
    }

    public SegmentedGenome getSegmentedGenome() { return segmentedGenome; }

    public AlleleFractionStateError error(final AlleleFractionState state) {
        final double averageMinorFractionError = IntStream.range(0, numSegments)
                .mapToDouble(s -> Math.abs(trueState.segmentMinorFraction(s) - state.segmentMinorFraction(s)))
                .average().getAsDouble();
        return new AlleleFractionStateError(averageMinorFractionError, trueState.meanBias() - state.meanBias(),
                trueState.biasVariance() - state.biasVariance(), trueState.outlierProbability() - state.outlierProbability());
    }

    public static final class AlleleFractionStateError {
        public final double averageMinorFractionError;
        public final double biasMeanError;
        public final double biasVarianceError;
        public final double outlierProbabilityError;

        public AlleleFractionStateError(final double averageMinorFractionError, final double biasMeanError,
                                        final double biasVarianceError, final double outlierProbabilityError) {
            this.averageMinorFractionError = averageMinorFractionError;
            this.biasMeanError = biasMeanError;
            this.biasVarianceError = biasVarianceError;
            this.outlierProbabilityError = outlierProbabilityError;
        }
    }

    public static final class PhaseIndicators extends ArrayList<AlleleFractionIndicator> {
        private static final long serialVersionUID = 60652L;
        public PhaseIndicators(final List<AlleleFractionIndicator> outlierIndicators) {
            super(new ArrayList<>(outlierIndicators));
        }
    }
}
