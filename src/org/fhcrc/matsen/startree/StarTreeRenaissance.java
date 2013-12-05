package org.fhcrc.matsen.startree;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import com.google.common.base.Function;
import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.*;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.app.beagle.evomodel.utilities.DnDsLogger;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.util.TaxonList;
import dr.evomodel.MSSD.ExponentialBranchLengthTreePrior;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.Logger;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.*;
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

    public static final int CHAIN_LENGTH = 20000;
    public static final int N_SAMPLES = 1000;
    public static final int SAMPLE_FREQ = CHAIN_LENGTH / N_SAMPLES;

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
        Preconditions.checkArgument(subsModels.size() == 3,
                "invalid number of substitution models: %d", subsModels.size());
        Preconditions.checkArgument(siteModels.size() == 3,
                "Invalid number of site models: %d", siteModels.size());

        // Patterns
        int maxIndex = alignment.getPatternCount();
        maxIndex -= maxIndex % 3;
        SitePatterns[] p = new SitePatterns[3];
        for (int i = 0; i < 3; i++) {
            p[i] = new SitePatterns(alignment, alignment, i, maxIndex - 1, 3, false, false);
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

        // DnDs Counts
        CodonPartitionedRobustCounting[] robustCounts = getCodonPartitionedRobustCountings(treeModel, treeLikelihoods);

        // Priors
        List<Likelihood> priors = new ArrayList<Likelihood>(1);
        DistributionLikelihood rootHeightPrior = new DistributionLikelihood(new ExponentialDistribution(10.0));
        rootHeightPrior.addData(treeModel.getRootHeightParameter());
        priors.add(rootHeightPrior);

        // Likelihood
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
        mcmc.init(options, like, operatorSchedule, new Logger[]{logger});
        mcmc.run();

        return twoTaxonResultOfTraces(formatter.getTraces(), p[0].getPatternCount());
    }

    private static TwoTaxonResult twoTaxonResultOfTraces(final List<Trace> traces, final int nCodons) {
        final int traceLength = traces.get(0).getValuesSize();
        final DoubleMatrix2D cn = new DenseDoubleMatrix2D(traceLength, nCodons), cs = new DenseDoubleMatrix2D(traceLength, nCodons),
                un = new DenseDoubleMatrix2D(traceLength, nCodons), us = new DenseDoubleMatrix2D(traceLength, nCodons);

        java.util.regex.Pattern p = java.util.regex.Pattern.compile("([CU])([NS])\\[(\\d+)\\]$");
        for (int i = 0; i < traceLength; i++) {
            for (final Trace trace : traces) {
                final String name = trace.getName();
                if("state".equals(name))
                    continue;
                final Matcher m = p.matcher(name);
                Preconditions.checkState(m.matches(), "%s does not match", name);

                final int pos = Integer.parseInt(m.group(3)) - 1;

                boolean isConditioned = name.substring(0, 1).equals(CodonPartitionedRobustCounting.SITE_SPECIFIC_PREFIX);
                boolean isNonsynonymous = name.substring(2, 3).equals(CodonLabeling.NON_SYN.getText());
                DoubleMatrix2D target;
                if (isConditioned && isNonsynonymous)
                    target = cn;
                else if (isConditioned)
                    target = cs;
                else if (isNonsynonymous)
                    target = un;
                else
                    target = us;

                @SuppressWarnings("unchecked")
                Trace<Double> dTrace = (Trace<Double>) trace;

                target.set(i, pos, dTrace.getValue(i));
            }
        }

        return new TwoTaxonResult(cn, cs, un, us);
    }

    private static CodonPartitionedRobustCounting[] getCodonPartitionedRobustCountings(final TreeModel treeModel, final AncestralStateBeagleTreeLikelihood[] treeLikelihoods) {
        CodonPartitionedRobustCounting[] robustCounts = new CodonPartitionedRobustCounting[2];

        String[] type = new String[]{"N", "S"};
        StratifiedTraitOutputFormat branchFormat = StratifiedTraitOutputFormat.SUM_OVER_SITES;
        StratifiedTraitOutputFormat logFormat = StratifiedTraitOutputFormat.SUM_OVER_SITES;
        for (int i = 0; i < 2; i++) {
            final String label = type[i];
            robustCounts[i] = new CodonPartitionedRobustCounting(label,
                    treeModel,
                    treeLikelihoods,
                    Codons.UNIVERSAL,
                    CodonLabeling.parseFromString(type[i]),
                    true,  // uniformization
                    true,  // external branches
                    true,  // internal branches
                    false, // unconditional per branch
                    false, // complete history
                    branchFormat,
                    logFormat);
        }
        return robustCounts;
    }

    private static AncestralStateBeagleTreeLikelihood[] getAncestralStateBeagleTreeLikelihoods(List<? extends SubstitutionModel> subsModels,
                                                                                               List<? extends SiteRateModel> siteModels,
                                                                                               SitePatterns[] p, StrictClockBranchRates branchRates, TreeModel treeModel) {
        AncestralStateBeagleTreeLikelihood[] treeLikelihoods = new AncestralStateBeagleTreeLikelihood[3];
        for (int i = 0; i < 3; i++) {
            treeLikelihoods[i] = new AncestralStateBeagleTreeLikelihood(p[i],
                    treeModel,
                    new HomogeneousBranchModel(subsModels.get(i), subsModels.get(i).getFrequencyModel()),
                    siteModels.get(i),
                    branchRates,
                    null,  // Tip states model
                    false, // Use ambiguities?
                    PartialsRescalingScheme.DEFAULT,
                    null,  // Partials restrictions
                    Nucleotides.INSTANCE,
                    String.format("CP%d.states", i + 1),
                    false, // Use MAP?
                    true); // Use ML?
        }
        return treeLikelihoods;
    }

    private static void createdNdSloggers(TreeModel treeModel, CodonPartitionedRobustCounting[] robustCounts, MCLogger logger) {
        final TreeTrait[] traits = new TreeTrait[4];
        int j = 0;
        for (final CodonPartitionedRobustCounting rc : robustCounts) {
            for (TreeTrait trait : rc.getTreeTraits()) {
                if(trait.getTraitName().matches("[CcUu]_[NnSs]")) {
                    traits[j++] = trait;
                }
            }
        }

        DnDsLogger dndsNLogger = new DnDsLogger("dndsN", treeModel, traits, false, false, true, false);
        logger.add(dndsNLogger);
        DnDsLogger dndsSLogger = new DnDsLogger("dndsS", treeModel, traits, false, false, true, true);
        logger.add(dndsSLogger);
    }

    public static void main(String... args) throws Exception {
        SAMFileReader reader = new SAMFileReader(new File("test.bam"));
        File fasta = new File("ighvdj.fasta");

        BufferedReader jsonReader = new BufferedReader(new FileReader("test.json"));
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
