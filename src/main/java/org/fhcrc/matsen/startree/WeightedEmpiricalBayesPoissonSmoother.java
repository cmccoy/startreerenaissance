package org.fhcrc.matsen.startree;

import com.google.common.base.Preconditions;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.moment.Variance;

/**
 * Created by cmccoy on 12/16/13.
 */
class WeightedEmpiricalBayesPoissonSmoother {
    private WeightedEmpiricalBayesPoissonSmoother() {
        Preconditions.checkArgument(false, "Do not construct.");
    }


    /**
     * This is equivalent @link dr.math.EmpiricalBayesPoissonSmoother, but sample mean and variance calculations are weighted
     *
     * For details, see:
     *
     * Lemey, Philippe, et al. "A counting renaissance: combining stochastic mapping and empirical Bayes to quickly detect amino acid sites under positive selection."
     * Bioinformatics 28.24 (2012): 3248-3256.
     *
     * @param values  Values to smooth
     * @param weights Weight of each value
     * @return Smoothed version
     */
    public static double[] smooth(final double[] values, final double[] weights) {
        Preconditions.checkNotNull(values, "Missing values");
        Preconditions.checkNotNull(weights, "Missing weights");
        Preconditions.checkArgument(values.length == weights.length,
                "value / weight lengths differ: %d vs %d",
                values.length, weights.length);

        // Transform so that the highest weight site has w=1.0
        final double maxWeight = StatUtils.max(weights);
        final double[] normWeights = new double[weights.length];
        for(int i = 0; i < weights.length; i++)
            normWeights[i] = weights[i] / maxWeight;

        // Lemey et. al. 2012 Equation 4
        final double mu = new Mean().evaluate(values, normWeights);
        final double sigma_sq = new Variance().evaluate(values, normWeights, mu);

        // Match moments of gamma poisson
        // Lemey et. al. 2012 Equation 5
        final double alpha = Math.pow(mu, 2) / (sigma_sq - mu),
                beta = mu / (sigma_sq - mu);

        final double[] result = new double[values.length];

        // Lemey et. al. 2012 Equation 7:
        if (sigma_sq > mu) {
            for (int i = 0; i < result.length; i++) {
                result[i] = (values[i] + alpha) / (1 + beta);
            }
        } else {
            java.util.Arrays.fill(result, mu);
        }

        return result;
    }
}
