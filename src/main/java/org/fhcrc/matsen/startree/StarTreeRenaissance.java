package org.fhcrc.matsen.startree;

import com.google.common.base.Preconditions;
import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.*;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.evomodel.utilities.DnDsLogger;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.util.Date;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.Logger;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.NegativeStatistic;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.RandomWalkOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.trace.Trace;
import dr.math.distributions.ExponentialDistribution;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.BlockRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.StatUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.locks.ReentrantLock;
import java.util.logging.Level;
import java.util.regex.Matcher;


/**
 * Created with IntelliJ IDEA.
 * User: cmccoy
 * Date: 12/2/13
 * Time: 12:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class StarTreeRenaissance {

    // Lock beagle instance creation-related actions - getting using instances that have been created from multiple threads
    // should present no problem.
    // TODO: what about releasing beagle instance?
    private static ReentrantLock beagleLock = new ReentrantLock();
    public static final int CHAIN_LENGTH = 20000;
    public static final int N_SAMPLES = 500;
    public static final int SAMPLE_FREQ = CHAIN_LENGTH / N_SAMPLES;

    private static final java.util.logging.Logger log = java.util.logging.Logger.getLogger("org.fhcrc.matsen.startree");

    private static final String ROOT_HEIGHT_NAME = "rootHeight";

    /**
     * Find start of first non-gap codon
     */
    private static int findOffset(final Alignment alignment) {
        Preconditions.checkState(alignment.getSequenceCount() == 2,
                "Unexpected number of sequences");
        final String query = alignment.getAlignedSequenceString(1);
        int offset = 0;
        for (int i = 0; i < query.length(); i += 3) {
            for (int j = i; j < i + 3; j++) {
                if (query.charAt(j) != '-')
                    return offset;
            }
            offset = i;
        }
        return 0;
    }

    /**
     * Generate some MCMC samples of synonymous / nonsynonymous substitutions both conditioned and unconditioned on data.
     * <p/>
     * Creates @link N_SAMPLES from a chain length of @link CHAIN_LENGTH
     *
     * @param alignment  An alignment with two taxa - no tree moves are performed!
     * @param subsModels One substitution model for each site.
     * @param siteModels One site model for each site.
     * @return A StarTreeTraces, with counts by iteration.
     * @throws Tree.MissingTaxonException
     */
    public static StarTreeTraces calculate(final Alignment alignment,
                                           final List<? extends SubstitutionModel> subsModels,
                                           final List<? extends SiteRateModel> siteModels) throws Tree.MissingTaxonException {
        return calculate(alignment, subsModels, siteModels, CHAIN_LENGTH, SAMPLE_FREQ);
    }

    private static double[] getCodonCoverage(final Alignment alignment) {
        Preconditions.checkNotNull(alignment);
        Preconditions.checkArgument(alignment.getSequenceCount() == 2,
                "Expected 2 sequences, got %s", alignment.getSequenceCount());

        final double[] result = new double[alignment.getPatternCount() / 3];
        final String qry = alignment.getAlignedSequenceString(1);

        for (int i = 0; i < result.length; i++) {
            result[i] = 0.0;
            for (int j = 0; j < 3; j++) {
                final char c = qry.charAt(3 * i + j);
                if (c != 'N' && c != '-') {
                    result[i] = 1;
                    continue;
                }
            }
        }

        Preconditions.checkState(StatUtils.max(result) == 1.0,
                "Unexpected maximum value: %f", StatUtils.max(result));
        return result;
    }

    /**
     * Generate some MCMC samples of synonymous / nonsynonymous substitutions both conditioned and unconditioned on data.
     *
     * @param alignment   An alignment with two taxa - no tree moves are performed!
     * @param subsModels  One substitution model for each site.
     * @param siteModels  One site model for each site.
     * @param chainLength Length of the MCMC chain
     * @param sampleEvery Sampling frequency
     * @return A StarTreeTraces, with counts by iteration.
     * @throws Tree.MissingTaxonException
     */
    public static StarTreeTraces calculate(final Alignment alignment,
                                           final List<? extends SubstitutionModel> subsModels,
                                           final List<? extends SiteRateModel> siteModels,
                                           final int chainLength,
                                           final int sampleEvery) throws Tree.MissingTaxonException {
        final long startTime = System.currentTimeMillis();
        Preconditions.checkArgument(subsModels.size() == 3,
                "invalid number of substitution models: %s", subsModels.size());
        Preconditions.checkArgument(siteModels.size() == 3,
                "Invalid number of site models: %s", siteModels.size());
        Preconditions.checkArgument(alignment.getSequenceCount() == 2,
                "Expected 2 sequences, got %s", alignment.getSequenceCount());

        // This is a hack, since this function often runs on a worker node, and we don't need to see citations for every pair
        final java.util.logging.Logger[] toDisable = new java.util.logging.Logger[]{
                java.util.logging.Logger.getLogger("dr.evomodel"),
                java.util.logging.Logger.getLogger("dr.app.beagle")
        };
        for (final java.util.logging.Logger l : toDisable) {
            l.setLevel(java.util.logging.Level.WARNING);
        }

        int minIndex = findOffset(alignment);
        log.log(Level.CONFIG, "offset: {0}", new Object[]{minIndex});
        Preconditions.checkState(minIndex % 3 == 0);

        // Patterns
        int maxIndex = alignment.getPatternCount();
        maxIndex -= maxIndex % 3;
        SitePatterns[] p = new SitePatterns[3];
        for (int i = 0; i < 3; i++) {
            p[i] = new SitePatterns(alignment, alignment, minIndex + i, maxIndex - 1, 3, false, false);
        }

        // Set tip dates
        alignment.getTaxon(0).setDate(new Date(0.00, Units.Type.SUBSTITUTIONS, false));
        alignment.getTaxon(1).setDate(new Date(0.00, Units.Type.SUBSTITUTIONS, false));

        // Coalescent model - this is just a shortcut to a starting tree.
        final ConstantPopulationModel coalModel = new ConstantPopulationModel(new Parameter.Default(1.0),
                dr.evolution.util.Units.Type.SUBSTITUTIONS);

        // Branch rates - fixed at 1.0 subs/site
        final StrictClockBranchRates branchRates = new StrictClockBranchRates(new Parameter.Default(1.0));

        // Tree model
        final CoalescentSimulator coalSim = new CoalescentSimulator();
        final TreeModel treeModel = new TreeModel("treeModel", coalSim.simulateTree(alignment, coalModel));

        final Parameter rootHeightParameter = treeModel.getRootHeightParameter();
        rootHeightParameter.setParameterValueQuietly(0, 1e-8);

        // We keep the reference sequence fixed at t = 0 (present),
        // and move the query sequence up and down in negative time.
        // This seems to work (in terms of the output trees), where I was unable to get co-varying root height and
        // reference height to work sensibly.
        final Parameter queryHeightParameter = treeModel.getLeafHeightParameter(treeModel.getExternalNode(1));
        queryHeightParameter.setId(ROOT_HEIGHT_NAME);
        queryHeightParameter.setParameterValueQuietly(0, -0.01);

        // Tree likelihoods
        AncestralStateBeagleTreeLikelihood[] treeLikelihoods = getAncestralStateBeagleTreeLikelihoods(subsModels, siteModels, p, branchRates, treeModel);

        // dNdS Counts
        RootConditionedCodonPartitionedRobustCounting[] robustCounts = getCodonPartitionedRobustCountings(treeModel, treeLikelihoods);

        // Priors
        List<Likelihood> priors = new ArrayList<Likelihood>(1);

        // The query height ends up being negative - so we just put an exponential prior on the negated value.
        DistributionLikelihood queryHeightPrior = new DistributionLikelihood(new ExponentialDistribution(10.0));
        queryHeightPrior.addData(new NegativeStatistic("negated." + queryHeightParameter.getParameterName(), queryHeightParameter));
        priors.add(queryHeightPrior);

        // Compound Likelihood
        List<Likelihood> likelihoods = new ArrayList<Likelihood>(4);
        likelihoods.add(new CompoundLikelihood(priors));
        likelihoods.get(0).setId("prior");
        likelihoods.addAll(Arrays.asList(treeLikelihoods));
        final CompoundLikelihood like = new CompoundLikelihood(likelihoods);
        like.setUsed();

        // Operators
        SimpleOperatorSchedule operatorSchedule = new SimpleOperatorSchedule();
        // Random walk on query height - this just moves it around in negative "time"
        operatorSchedule.addOperator(new RandomWalkOperator(
                queryHeightParameter,
                1.0,
                RandomWalkOperator.BoundaryCondition.reflecting,
                1.0,
                CoercionMode.DEFAULT));

        // log to a big array
        ArrayLogFormatter formatter = new ArrayLogFormatter(false);
        MCLogger logger = new MCLogger(formatter, sampleEvery, false);
        createdNdSloggers(treeModel, robustCounts, logger);
        logger.add(new NegativeStatistic(ROOT_HEIGHT_NAME, queryHeightParameter));

        // Actually run the MCMC
        MCMCOptions options = new MCMCOptions(chainLength);
        MCMC mcmc = new MCMC("mcmc");
        mcmc.setShowOperatorAnalysis(false);
        mcmc.init(options, like, operatorSchedule, new Logger[]{logger});
        mcmc.run();

        final StarTreeTraces result = twoTaxonResultOfTraces(formatter.getTraces(), p[0].getPatternCount(), minIndex / 3, getCodonCoverage(alignment));

        final long endTime = System.currentTimeMillis();

        log.log(Level.FINE, String.format("Sampled in in %d ms", endTime - startTime));

        return result;
    }

    /**
     * Pull the statistics we're interested out of a collection of traces
     */
    private static StarTreeTraces twoTaxonResultOfTraces(final List<Trace> traces, final int nCodons, final int offset, final double[] coverage) {
        final int traceLength = traces.get(0).getValuesSize();
        final RealMatrix cn = new BlockRealMatrix(traceLength, offset + nCodons),
                cs = new BlockRealMatrix(traceLength, offset + nCodons),
                un = new BlockRealMatrix(traceLength, offset + nCodons),
                us = new BlockRealMatrix(traceLength, offset + nCodons),
                bl = new BlockRealMatrix(traceLength, offset + nCodons);

        final RealVector state = new ArrayRealVector(traceLength);

        java.util.regex.Pattern p = java.util.regex.Pattern.compile("[CU][NS]\\[(\\d+)\\]$");
        for (int i = 0; i < traceLength; i++) {
            for (final Trace trace : traces) {
                final String name = trace.getName();
                if (name.equals("state")) {
                    @SuppressWarnings("unchecked")
                    final Trace<Double> dTrace = (Trace<Double>) trace;
                    state.setEntry(i, dTrace.getValue(i));
                    continue;
                } else if (name.equals(ROOT_HEIGHT_NAME)) {
                    @SuppressWarnings("unchecked")
                    final Trace<Double> dTrace = (Trace<Double>) trace;

                    final double d = dTrace.getValue(i);
                    for (int j = 0; j < bl.getColumnDimension(); j++) {
                        bl.setEntry(i, j, coverage[j] > 0 ? d : 0);
                    }
                    continue;
                }
                final Matcher m = p.matcher(name);
                Preconditions.checkState(m.matches(), "%s does not match", name);

                final int pos = Integer.parseInt(m.group(1)) - 1 + offset;

                final RealMatrix target;
                if (name.startsWith("CN"))
                    target = cn;
                else if (name.startsWith("CS"))
                    target = cs;
                else if (name.startsWith("UN"))
                    target = un;
                else if (name.startsWith("US"))
                    target = us;
                else
                    throw new IllegalStateException("unknown trace: " + name);

                @SuppressWarnings("unchecked")
                Trace<Double> dTrace = (Trace<Double>) trace;

                target.setEntry(i, pos, dTrace.getValue(i));
            }
        }

        return new StarTreeTraces(state, cn, cs, un, us, coverage, bl);
    }

    private static RootConditionedCodonPartitionedRobustCounting[] getCodonPartitionedRobustCountings(final TreeModel treeModel, final AncestralStateBeagleTreeLikelihood[] treeLikelihoods) {
        RootConditionedCodonPartitionedRobustCounting[] robustCounts = new RootConditionedCodonPartitionedRobustCounting[2];

        String[] type = new String[]{"S", "N"};
        StratifiedTraitOutputFormat branchFormat = StratifiedTraitOutputFormat.SUM_OVER_SITES;
        StratifiedTraitOutputFormat logFormat = StratifiedTraitOutputFormat.SUM_OVER_SITES;
        for (int i = 0; i < 2; i++) {
            robustCounts[i] = new RootConditionedCodonPartitionedRobustCounting(
                    type[i],
                    treeModel,
                    treeLikelihoods,
                    Codons.UNIVERSAL,
                    CodonLabeling.parseFromString(type[i]),
                    true,  // uniformization
                    true,  // external branches
                    true,  // internal branches
                    true,  // unconditional per branch
                    false, // complete history
                    branchFormat,
                    logFormat);
        }
        return robustCounts;
    }

    private static AncestralStateBeagleTreeLikelihood[] getAncestralStateBeagleTreeLikelihoods(List<? extends SubstitutionModel> subsModels,
                                                                                               List<? extends SiteRateModel> siteModels,
                                                                                               SitePatterns[] p,
                                                                                               StrictClockBranchRates branchRates,
                                                                                               TreeModel treeModel) {
        AncestralStateBeagleTreeLikelihood[] treeLikelihoods = new AncestralStateBeagleTreeLikelihood[3];
        try {
            beagleLock.lock();
            for (int i = 0; i < 3; i++) {
                treeLikelihoods[i] = new AncestralStateBeagleTreeLikelihood(p[i],
                        treeModel,
                        new HomogeneousBranchModel(subsModels.get(i), subsModels.get(i).getFrequencyModel()),
                        siteModels.get(i),
                        branchRates,
                        null,  // Tip states model
                        false, // Use ambiguities?
                        PartialsRescalingScheme.DELAYED,
                        null,  // Partials restrictions
                        Nucleotides.INSTANCE,
                        String.format("CP%d.states", i + 1),
                        false, // Use MAP?
                        true); // Use ML?
            }
        } finally {
            beagleLock.unlock();
        }
        return treeLikelihoods;
    }

    private static void createdNdSloggers(TreeModel treeModel, RootConditionedCodonPartitionedRobustCounting[] robustCounts, MCLogger logger) {
        // expected trait order: c_S u_S c_N u_N
        final TreeTrait[] traits = new TreeTrait[4];
        final String[] expectedNames = new String[]{"c_S", "u_S", "c_N", "u_N"};
        for (final RootConditionedCodonPartitionedRobustCounting rc : robustCounts) {
            for (TreeTrait trait : rc.getTreeTraits()) {
                if (trait.getTraitName().matches("[cu]_[NS]")) {
                    final String name = trait.getTraitName();
                    // Generate index to match expected order above
                    final int idx = (name.endsWith("S") ? 0 : 2) + (name.startsWith("c") ? 0 : 1);
                    Preconditions.checkState(traits[idx] == null);
                    traits[idx] = trait;
                }
            }
        }
        for (int i = 0; i < traits.length; i++) {
            Preconditions.checkState(traits[i] != null, "Missing trait %d", i);
            Preconditions.checkState(traits[i].getTraitName().equals(expectedNames[i]),
                    "Unexpected name: got %s, expected %s", traits[i].getTraitName(), expectedNames[i]);
            Preconditions.checkState(TreeTrait.Intent.WHOLE_TREE.equals(traits[i].getIntent()),
                    "Unsupported intent in trait %s: %s (need WHOLE_TREE)",
                    traits[i].getTraitName(),
                    traits[i].getIntent());
        }


        logger.add(new DnDsLogger("dndsN", treeModel, traits, false, false, true, false));
        logger.add(new DnDsLogger("dndsS", treeModel, traits, false, false, true, true));
    }
}
