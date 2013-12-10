package org.fhcrc.matsen.startree;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.RCDoubleMatrix2D;
import com.google.common.base.Function;
import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.CodonLabeling;
import dr.app.beagle.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.app.beagle.evomodel.substmodel.StratifiedTraitOutputFormat;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.evomodel.utilities.DnDsLogger;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
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
import dr.inference.model.Parameter;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.trace.Trace;
import dr.math.distributions.ExponentialDistribution;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.locks.ReentrantLock;
import java.util.logging.Level;
import java.util.regex.Matcher;


// TODOs
// * Run for each pair, aggregate matrices (spark?)
// * Verify that no smoothing is being performed

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
    public static final int N_SAMPLES = 1000;
    public static final int SAMPLE_FREQ = CHAIN_LENGTH / N_SAMPLES;

    private static final java.util.logging.Logger log = java.util.logging.Logger.getLogger("org.fhcrc.matsen.startree");

    /**
     * Generate some MCMC samples of synonymous / nonsynonymous substitutions both conditioned and unconditioned on data.
     *
     * @param alignment  An alignment with two taxa - no tree moves are performed!
     * @param subsModels One substitution model for each site.
     * @param siteModels One site model for each site.
     * @return A TwoTaxonResult, with counts by iteration.
     * @throws Tree.MissingTaxonException
     */
    public static TwoTaxonResult calculate(final Alignment alignment,
                                           final List<? extends SubstitutionModel> subsModels,
                                           final List<? extends SiteRateModel> siteModels) throws Tree.MissingTaxonException {
        return calculate(alignment, subsModels, siteModels, CHAIN_LENGTH, SAMPLE_FREQ);
    }

    /**
     * Find start of first non-gap codon
     */
    private static int findOffset(final Alignment alignment) {
      Preconditions.checkState(alignment.getSequenceCount() == 2,
          "Unexpected number of sequences");
      final String query = alignment.getAlignedSequenceString(1);
      int offset = 0;
      for(int i = 0; i < query.length(); i += 3) {
        for(int j = i; j < i + 3; j++) {
          if(query.charAt(j) != '-')
            return offset;
        }
        offset = i;
      }
      return 0;
    }

    /**
     * Generate some MCMC samples of synonymous / nonsynonymous substitutions both conditioned and unconditioned on data.
     *
     * @param alignment  An alignment with two taxa - no tree moves are performed!
     * @param subsModels One substitution model for each site.
     * @param siteModels One site model for each site.
     * @return A TwoTaxonResult, with counts by iteration.
     * @throws Tree.MissingTaxonException
     */
    public static TwoTaxonResult calculate(final Alignment alignment,
                                           final List<? extends SubstitutionModel> subsModels,
                                           final List<? extends SiteRateModel> siteModels,
                                           final int chainLength,
                                           final int sampleEvery) throws Tree.MissingTaxonException {
        final long startTime = System.currentTimeMillis();
        Preconditions.checkArgument(subsModels.size() == 3,
                "invalid number of substitution models: %d", subsModels.size());
        Preconditions.checkArgument(siteModels.size() == 3,
                "Invalid number of site models: %d", siteModels.size());

        // This is a hack, since this function often runs on a worker node, and we don't need to see citations for every pair
        java.util.logging.Logger.getLogger("dr.evomodel").setLevel(java.util.logging.Level.WARNING);
        java.util.logging.Logger.getLogger("dr.app.beagle").setLevel(java.util.logging.Level.WARNING);

        int minIndex = findOffset(alignment);
        log.log(Level.CONFIG, "offset: %d", new Object[] { minIndex });
        Preconditions.checkState(minIndex % 3 == 0);

        // Patterns
        int maxIndex = alignment.getPatternCount();
        maxIndex -= maxIndex % 3;
        SitePatterns[] p = new SitePatterns[3];
        for (int i = 0; i < 3; i++) {
            p[i] = new SitePatterns(alignment, alignment, minIndex + i, maxIndex - 1, 3, false, false);
        }

        // Coalescent model
        final ConstantPopulationModel coalModel = new ConstantPopulationModel(new Parameter.Default(1.0),
                dr.evolution.util.Units.Type.SUBSTITUTIONS);

        // Branch rates - fixed at 1.0 subs/site
        final StrictClockBranchRates branchRates = new StrictClockBranchRates(new Parameter.Default(1.0));

        // Tree model
        CoalescentSimulator coalSim = new CoalescentSimulator();
        final TreeModel treeModel = new TreeModel("treeModel", coalSim.simulateTree(alignment, coalModel));
        treeModel.getRootHeightParameter().setParameterValueQuietly(0, 0.05);

        // Tree likelihoods
        AncestralStateBeagleTreeLikelihood[] treeLikelihoods = getAncestralStateBeagleTreeLikelihoods(subsModels, siteModels, p, branchRates, treeModel);

        // dNdS Counts
        CodonPartitionedRobustCounting[] robustCounts = getCodonPartitionedRobustCountings(treeModel, treeLikelihoods);

        // Priors
        List<Likelihood> priors = new ArrayList<Likelihood>(1);
        DistributionLikelihood rootHeightPrior = new DistributionLikelihood(new ExponentialDistribution(10.0));
        rootHeightPrior.addData(treeModel.getRootHeightParameter());
        priors.add(rootHeightPrior);

        // Compound Likelihood
        List<Likelihood> likelihoods = new ArrayList<Likelihood>(3);
        likelihoods.add(new CompoundLikelihood(priors));
        likelihoods.addAll(Arrays.asList(treeLikelihoods));
        final CompoundLikelihood like = new CompoundLikelihood(likelihoods);
        like.setUsed();

        // Operators
        SimpleOperatorSchedule operatorSchedule = new SimpleOperatorSchedule();
        operatorSchedule.addOperator(new ScaleOperator(treeModel.getRootHeightParameter(), 0.75));

        // log
        ArrayLogFormatter formatter = new ArrayLogFormatter(false);
        MCLogger logger = new MCLogger(formatter, sampleEvery, false);
        createdNdSloggers(treeModel, robustCounts, logger);

        MCMCOptions options = new MCMCOptions(chainLength);
        MCMC mcmc = new MCMC("mcmc");
        mcmc.setShowOperatorAnalysis(false);
        mcmc.init(options, like, operatorSchedule, new Logger[]{logger});
        mcmc.run();

        final TwoTaxonResult result = twoTaxonResultOfTraces(formatter.getTraces(), p[0].getPatternCount(), minIndex / 3);

        final long endTime = System.currentTimeMillis();

        log.log(Level.INFO, String.format("Sampled in in %d ms", endTime - startTime));

        return result;
    }

    /**
     * Pull the statistics we're interested out of a collection of traces
     */
    private static TwoTaxonResult twoTaxonResultOfTraces(final List<Trace> traces, final int nCodons, final int offset) {
        final int traceLength = traces.get(0).getValuesSize();
        final DoubleMatrix2D cn = new RCDoubleMatrix2D(traceLength, offset + nCodons),
                cs = new RCDoubleMatrix2D(traceLength, offset + nCodons),
                un = new RCDoubleMatrix2D(traceLength, offset + nCodons),
                us = new RCDoubleMatrix2D(traceLength, offset + nCodons);

        final DoubleMatrix1D totalN = new DenseDoubleMatrix1D(traceLength);
        final DoubleMatrix1D totalS = new DenseDoubleMatrix1D(traceLength);
        final DoubleMatrix1D state = new DenseDoubleMatrix1D(traceLength);

        java.util.regex.Pattern p = java.util.regex.Pattern.compile("([CU])([NS])\\[(\\d+)\\]$");
        for (int i = 0; i < traceLength; i++) {
            for (final Trace trace : traces) {
                final String name = trace.getName();
                if(name.matches("state|total_[NS]")) {
                    final DoubleMatrix1D target;
                    if("state".equals(name))
                        target = state;
                    else if("total_N".equals(name))
                        target = totalN;
                    else if("total_S".equals(name))
                        target = totalS;
                    else
                        throw new IllegalStateException("unknown trace: " + name);
                    @SuppressWarnings("unchecked")
                    final Trace<Double> dTrace = (Trace<Double>) trace;
                    target.setQuick(i, dTrace.getValue(i));
                    continue;
                }
                final Matcher m = p.matcher(name);
                Preconditions.checkState(m.matches(), "%s does not match", name);

                final int pos = Integer.parseInt(m.group(3)) - 1 + offset;

                final DoubleMatrix2D target;
                if(name.startsWith("CN"))
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

                target.set(i, pos, dTrace.getValue(i));
            }
        }

        final DoubleMatrix1D cnSum = new DenseDoubleMatrix1D(traceLength),
                             csSum = new DenseDoubleMatrix1D(traceLength),
                             unSum = new DenseDoubleMatrix1D(traceLength),
                             usSum = new DenseDoubleMatrix1D(traceLength);
        for(int i = 0; i < cn.columns(); i++) {
          cnSum.assign(cn.viewColumn(i), cern.jet.math.Functions.plus);
          csSum.assign(cs.viewColumn(i), cern.jet.math.Functions.plus);
          unSum.assign(un.viewColumn(i), cern.jet.math.Functions.plus);
          usSum.assign(us.viewColumn(i), cern.jet.math.Functions.plus);
        }

        return new TwoTaxonResult(state, cn, cs, un, us, totalN, totalS);
    }

    private static CodonPartitionedRobustCounting[] getCodonPartitionedRobustCountings(final TreeModel treeModel, final AncestralStateBeagleTreeLikelihood[] treeLikelihoods) {
        CodonPartitionedRobustCounting[] robustCounts = new CodonPartitionedRobustCounting[2];

        String[] type = new String[]{"S", "N"};
        StratifiedTraitOutputFormat branchFormat = StratifiedTraitOutputFormat.SUM_OVER_SITES;
        StratifiedTraitOutputFormat logFormat = StratifiedTraitOutputFormat.SUM_OVER_SITES;
        for (int i = 0; i < 2; i++) {
            robustCounts[i] = new CodonPartitionedRobustCounting(
                    type[i],
                    treeModel,
                    treeLikelihoods,
                    Codons.UNIVERSAL,
                    CodonLabeling.parseFromString(type[i]),
                    true,  // uniformization
                    true,  // external branches
                    true,  // internal branches
                    true, // unconditional per branch
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

    private static void createdNdSloggers(TreeModel treeModel, CodonPartitionedRobustCounting[] robustCounts, MCLogger logger) {
        final TreeTrait[] traits = new TreeTrait[4];
        int j = 0;
        for (final CodonPartitionedRobustCounting rc : robustCounts) {
            logger.add(rc);
            for (TreeTrait trait : rc.getTreeTraits()) {
                if(trait.getTraitName().matches("[CcUu]_[NnSs]")) {
                    traits[j++] = trait;
                }
            }
        }

        logger.add(new DnDsLogger("dndsN", treeModel, traits, false, false, true, false));
        logger.add(new DnDsLogger("dndsS", treeModel, traits, false, false, true, true));
    }

    public static void main(String... args) throws Exception {
        SAMFileReader reader = new SAMFileReader(new File("simulate/merged.bam"));
        File fasta = new File("simulate/merged.fasta");

        BufferedReader jsonReader = new BufferedReader(new FileReader("testdata/test.json"));
        List<HKYModelParser.HKYAndRate> mRates = HKYModelParser.substitutionModel(jsonReader);

        final List<SubstitutionModel> models = new ArrayList<SubstitutionModel>();
        final List<SiteRateModel> rates = new ArrayList<SiteRateModel>();
        for (int i = 0; i < mRates.size(); i++) {
            models.add(mRates.get(i).getModel());
            GammaSiteRateModel r = new GammaSiteRateModel(String.format("rate%d", i));
            r.setMu(mRates.get(i).getRate());
            rates.add(r);
        }

        final Map<String, byte[]> references = SAMUtils.readAllFasta(fasta);

        Iterable<Alignment> alignments = Iterables.transform(reader, new Function<SAMRecord, Alignment>() {
            @Override
            public Alignment apply(SAMRecord samRecord) {
                final byte[] ref = references.get(samRecord.getReferenceName());
                Preconditions.checkNotNull(ref, "no reference for %s", samRecord.getReferenceName());
                return SAMBEASTUtils.alignmentOfRecord(samRecord, ref);
            }
        });

        Iterable<TwoTaxonResult> results = Iterables.transform(alignments, new Function<Alignment, TwoTaxonResult>() {
            @Override
            public TwoTaxonResult apply(Alignment a) {
                try {
                    return calculate(a, models, rates);
                } catch (Tree.MissingTaxonException exception) {
                    throw new RuntimeException(exception);
                }
            }
        });
    }
}
