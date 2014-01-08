package org.fhcrc.matsen.startree;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import dr.math.Poisson;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * <p>
 * Alternative version of @link dr.math.EmpiricalBayesPoissonSmoother, where coverage across positions is presumed to vary.
 * <p/>
 * <p>
 * we assume that <c>C_\ell \sim \mathrm{Poisson}(\lambda_\ell t_\ell)</c>,
 * where <c>t_\ell = \sum_i t_{\ell i}</c>, the sum of the branch lengths of sequences with a codon present at position <c>\ell</c>.
 * </p>
 * Created by cmccoy on 12/27/13.
 */
public class TLambdaPoissonSmoother {

    private static final Logger logger = Logger.getLogger(TLambdaPoissonSmoother.class.getName());

    private final static double BRENT_TOL = 1e-3;
    private final static int ADDL_INTERPOLATION_PTS = 2;

    private TLambdaPoissonSmoother() {
        Preconditions.checkArgument(false, "Do not construct.");
    }

    /**
     * Smooth counts <c>c</c> <i>without</i> sampling.
     *
     * @param c Substitution counts
     * @param t Total branch lengths for each entry in <c>c</c>
     * @return Smoothed rates
     */
    public static double[] smooth(final double[] c, final double[] t) {
        return smooth(c, t, false);
    }

    private static double[] smoothPoissonOnly(final double[] c, final double[] t, final boolean sample) {
        final double[] smoothed = new double[c.length];

        final double lambda = fitLambda(c, t);

        for (int i = 0; i < c.length; i++) {
            if (sample) {
                smoothed[i] = Poisson.nextPoisson(lambda * t[i]);
            } else {
                smoothed[i] = lambda;
            }
        }

        return smoothed;
    }

    /**
     * Smooth the counts <c>c</c>, producing per-site rates
     * <p/>
     * In the event that the variance of lambda is less than the mean, we fit:
     * <p/>
     * <code>
     * c_i \sim Poisson(t_i * lambda)
     * </code>
     * <p/>
     * Where lambda is estimated from the data
     *
     * @param c      Substitution counts
     * @param t      Total branch lengths for each entry in <c>c</c>
     * @param sample Should rates be sampled from the posterior distribution? If not, posterior mean is used.
     * @return Smoothed rates
     */
    public static double[] smooth(final double[] c, final double[] t, final boolean sample) {
        Preconditions.checkNotNull(c);
        Preconditions.checkNotNull(t);
        Preconditions.checkArgument(c.length == t.length,
                "Non-matching array lengths: %s vs %s", c.length, t.length);

        // Empirical Bayes
        logger.log(Level.FINE, "Starting Empirical Bayes estimation of alpha, beta.");

        final double[] smoothed = new double[c.length];

        final PointValuePair result;
        try {
            result = estimateAlphaBeta(c, t);
        } catch (org.apache.commons.math3.exception.MathIllegalStateException e) {
            final double mean = new Mean().evaluate(c, t);
            final double var = new Variance().evaluate(c, t, mean);
            logger.log(Level.WARNING, "Optimization failed. mean: {0} var: {1}. Using Poisson.\n{2}",
                    new Object[]{mean, var, e});

            return smoothPoissonOnly(c, t, sample);
        }

        final double alpha = result.getPoint()[0],
                beta = result.getPoint()[1];
        logger.log(Level.FINE, "MLE: alpha={0} beta={1} logL={2}",
                new Object[]{alpha, beta, result.getValue()});

        for (int i = 0; i < c.length; i++) {
            // posterior mean = alpha / beta
            if (sample) {
                final double shape = c[i] + alpha;
                final double scale = 1 / (t[i] + beta);
                smoothed[i] = dr.math.distributions.GammaDistribution.nextGamma(shape, scale);
            } else {
                smoothed[i] = (c[i] + alpha) / (t[i] + beta);
            }
        }

        return smoothed;
    }

    /**
     * Fits:
     * <code>
     * c_i \sim Poisson(t_i * lambda)
     * </code>
     * <p/>
     * Using Brent's method.
     *
     * @param c counts
     * @param t branch lengths
     * @return estimate of lambda
     */
    @VisibleForTesting
    static double fitLambda(final double[] c, final double[] t) {
        final BrentOptimizer optimizer = new BrentOptimizer(BRENT_TOL, BRENT_TOL);

        final UnivariateFunction fn = new UnivariateFunction() {
            @Override
            public double value(double x) {
                // poisson pmf = lambda^k / k! * e^-lambda
                //             = lambda^k / Gamma(k + 1) * e^-lambda
                // In our case, lambda = x * t[i]
                // log likelihood:
                // log(x * ti) * ci - logGamma(ci+1) - x * ti
                double result = 0;
                for (int i = 0; i < c.length; i++) {
                    final double lambda = Math.max(x * t[i], 1e-6);
                    PoissonDistribution p = new PoissonDistribution(lambda);
                    result += Math.log(p.probability((int) c[i]));
                }

                return result;
            }
        };

        final UnivariatePointValuePair result = optimizer.optimize(
                new MaxEval(100),
                new UnivariateObjectiveFunction(fn),
                GoalType.MAXIMIZE,
                new SearchInterval(0.001, 1000));

        logger.log(Level.FINE, "Brent finished in {0} steps",
                optimizer.getIterations());

        return result.getPoint();
    }

    /**
     * Estimate alpha, beta given c, t
     *
     * @param c Counts
     * @param t Branch lengths
     * @return Point, containing (alpha, beta), and associated log-likelihood
     */
    @VisibleForTesting
    static PointValuePair estimateAlphaBeta(final double[] c, final double[] t) {
        // Optimize. point[0] = alpha, point[1] = beta
        final int dim = 2;
        BOBYQAOptimizer optimizer = new BOBYQAOptimizer(2 * dim + ADDL_INTERPOLATION_PTS);

        final MultivariateFunction fn = new MultivariateFunction() {
            @Override
            public double value(double[] point) {
                Preconditions.checkArgument(point.length == 2,
                        "Invalid data size: %s", point.length);
                return logLikelihood(point[0], point[1], c, t);
            }
        };

        // Start with method-of-moments
        final double mean = StatUtils.mean(c),
                var = StatUtils.variance(c);

        double mom_alpha, mom_beta;
        if (var > mean) {
            mom_alpha = (mean * mean) / (var - mean);
            mom_beta = mean / (var - mean);
        } else {
            mom_alpha = 1;
            mom_beta = 1;
        }

        // Bounds - both must be positive
        final double lowerBounds[] = new double[]{1e-7, 1e-7};
        final double upperBounds[] = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        final double initial[] = new double[]{mom_alpha, mom_beta};
        try {
            PointValuePair result = optimizer.optimize(
                    new MaxEval(2000),
                    new InitialGuess(initial),
                    GoalType.MAXIMIZE,
                    new ObjectiveFunction(fn),
                    new SimpleBounds(lowerBounds, upperBounds));
            logger.log(Level.FINE, "BOBYQA finished in {0} iterations", optimizer.getIterations());
            return result;
        } catch (Exception e) {
//            logger.log(Level.SEVERE,
//                    "Optimization failed [mom_alpha: {0}, mom_beta: {1}]\nc={2}\nt={3}",
//                    new Object[]{mom_alpha, mom_beta,
//                            java.util.Arrays.toString(c),
//                            java.util.Arrays.toString(t)});
            throw e;
        }
    }

    /**
     * Returns the log likelihood of (c | t, alpha, beta)
     * <p/>
     * <c>\mathcal L = \left(\frac{\beta^\alpha}{\Gamma(\alpha)}\right)^L \prod_\ell \frac{t_\ell^{C_\ell}}{\Gamma(C_\ell+1)} \frac{\Gamma(C_\ell + \alpha)}{(t_\ell + \beta)^{C_\ell + \alpha}}</c>
     *
     * @param alpha Alpha parameter of gamma distribution
     * @param beta  Beta parameter of gamma distribution
     * @param c     Counts
     * @param t     Branch lengths
     * @return log-likelihood
     */
    @VisibleForTesting
    static double logLikelihood(final double alpha, final double beta,
                                final double[] c, final double[] t) {
        Preconditions.checkNotNull(c);
        Preconditions.checkNotNull(t);
        Preconditions.checkArgument(c.length == t.length,
                "Non-matching array lengths: %s vs %s", c.length, t.length);

        double result = 0.0;
        final int n = c.length;

        // \mathcal L = \left(\frac{\beta^\alpha}{\Gamma(\alpha)}\right)^L \prod_\ell \frac{t_\ell^{C_\ell}}{\Gamma(C_\ell+1)} \frac{\Gamma(C_\ell + \alpha)}{(t_\ell + \beta)^{C_\ell + \alpha}} $$

        result = 0;

        // Over sites
        for (int i = 0; i < n; i++) {
            final double ci = c[i], ti = t[i];

            // Special case for zero branch length (no coverage)
            if (ti <= 1e-8) {
                continue;
            }

            // Per Wolfram Alpha,
            // log(b^a×(t^c/(Gamma(c+1)))/(Gamma(a))×(Gamma(c+a))/(t+b)^(c+a)) simplifies to
            // (-a-c) log(b+t)+a log(b)+log(Gamma(a+c))-log(Gamma(a))+c log(t)-log(Gamma(c+1))
            double x = (-alpha - ci) * Math.log(beta + ti) +
                    alpha * Math.log(beta) +
                    Gamma.logGamma(alpha + ci) -
                    Gamma.logGamma(alpha) +
                    ci * Math.log(ti) -
                    Gamma.logGamma(ci + 1);
            result += x;
        }

        logger.log(Level.FINE, "alpha={0} beta={1} LL={2}",
                new Object[]{alpha, beta, result});

        return result;
    }
}
