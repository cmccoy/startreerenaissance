package org.fhcrc.matsen.startree;

import com.google.common.base.Preconditions;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Alternative version of @link dr.math.EmpiricalBayesPoissonSmoother, where coverage across positions is presumed to vary.
 *
 * we assume that <c>C_\ell \sim \mathrm{Poisson}(\lambda_\ell t_\ell)</c>,
 * where <c>t_\ell = \sum_i t_{\ell i}</c>, the sum of the branch lengths of sequences with a codon present at position <c>\ell</c>.
 *
 * Created by cmccoy on 12/27/13.
 */
public class TLambdaPoissonSmoother {

    private static final Logger logger = Logger.getLogger(TLambdaPoissonSmoother.class.getName());

    private final static int ADDL_INTERPOLATION_PTS = 1;

    /**
     * Smooth the counts <c>c</c>,
     * @param c
     * @param t
     * @return
     */
    public static double[] smooth(final double[] c, final double[] t) {
        Preconditions.checkNotNull(c);
        Preconditions.checkNotNull(t);
        Preconditions.checkArgument(c.length == t.length,
                "Non-matching array lengths: %s vs %s", c.length, t.length);

        // Empirical Bayes
        logger.info("Starting Empirical Bayes estimation of alpha, beta.");
        PointValuePair result = estimateAlphaBeta(c, t);
        final double alpha = result.getPoint()[0],
                beta = result.getPoint()[1];
        logger.log(Level.INFO, "ML result: alpha={0} beta={1} ll={2}",
                new Object[]{ alpha, beta, result.getValue()});

        final double[] smoothed = new double[c.length];

        // TODO: support random draws
        for(int i = 0; i < c.length; i++) {
            // posterior mean = alpha /
            smoothed[i] = (c[i] + alpha) / (t[i] + beta);
        }

        return smoothed;
    }

    private static PointValuePair estimateAlphaBeta(final double[] c, final double[] t) {
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

        final double lowerBounds[] = new double[] {0.0, 0.0};
        final double upperBounds[] = new double[] {5000, 5000};
        PointValuePair result = optimizer.optimize(
                new MaxEval(1000),
                GoalType.MAXIMIZE,
                new ObjectiveFunction(fn),
                new SimpleBounds(lowerBounds, upperBounds),
                new InitialGuess(new double[]{1, 1}));
        return result;
    }

    /**
     * Returns the log likelihood of (c | t, alpha, beta)
     *
     * <c>\mathcal L = \left(\frac{\beta^\alpha}{\Gamma(\alpha)}\right)^L \prod_\ell \frac{t_\ell^{C_\ell}}{\Gamma(C_\ell+1)} \frac{\Gamma(C_\ell + \alpha)}{(t_\ell + \beta)^{C_\ell + \alpha}}</c>
     *
     * @param alpha Alpha parameter of gamma distribution
     * @param beta Beta parameter of gamma distribution
     * @param c Counts
     * @param t Branch lengths
     * @return log-likelihood
     */
    private static double logLikelihood(final double alpha, final double beta,
                                        final double[] c, final double[] t) {
        Preconditions.checkNotNull(c);
        Preconditions.checkNotNull(t);
        Preconditions.checkArgument(c.length == t.length,
                "Non-matching array lengths: %s vs %s", c.length, t.length);

        double result = 0.0;
        final int n = c.length;

        // \mathcal L = \left(\frac{\beta^\alpha}{\Gamma(\alpha)}\right)^L \prod_\ell \frac{t_\ell^{C_\ell}}{\Gamma(C_\ell+1)} \frac{\Gamma(C_\ell + \alpha)}{(t_\ell + \beta)^{C_\ell + \alpha}} $$



        // log(((alpha ^ beta / Gamma(alpha))^N)
        // (log(alpha) * beta - log(Gamma(alpha)) * N
        result = (Math.log(alpha) * beta - Gamma.logGamma(alpha)) * n;

        // Over sites
        for(int i = 0; i < n; i++) {
            // t_l ^ C_l / Gamma(C_l + 1)
            double x = Math.log(t[i]) * c[i] - Gamma.logGamma(c[i] + 1.0);

            // Gamma(C_l + alpha) / (t_l + beta) ^ (C_l + alpha)
            x += Gamma.logGamma(c[i] + alpha);
            x -= Math.log(t[i] + beta) * (c[i] + alpha);

            result += x;
        }

        logger.log(Level.INFO, "LL = {0}", result);
        return result;
    }
}
