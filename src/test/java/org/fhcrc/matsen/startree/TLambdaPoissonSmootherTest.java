package org.fhcrc.matsen.startree;

import dr.math.EmpiricalBayesPoissonSmoother;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.StatUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.FileWriter;
import java.util.Arrays;

/**
 * Created by cmccoy on 12/27/13.
 *
 * <code>
 * set.seed(1)
 *
 * n <- 1000
 *
 * alpha <- 1.0
 * beta <- 2.0
 * rg <- rgamma(n, alpha, rate=beta)
 * tl <- rep(1.0, n)
 * rp <- unlist(lapply(seq_along(tl), function(i) rpois(1, tl[i] * rg[i])))
 *
 * cat(paste(rp, collapse=', '), '\n')
 * </code>
 */
public class TLambdaPoissonSmootherTest {
    public final double TOL = 3e-2;

    @Test
    public void testSmooth() throws Exception {
        final double[] arr = new double[] {0, 1, 3, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 4, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 2, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 1, 1, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 2, 0, 3, 1, 0, 0, 0, 0, 2, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 3, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 1, 1, 1, 0, 0, 2, 2, 3, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 6, 1, 0, 3, 1, 1, 3, 0, 2, 0, 1, 0, 0, 0, 2, 1, 4, 1, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 0, 0, 1, 1, 5, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 2, 1, 0, 2, 1, 0, 0, 0, 1, 0, 0, 2, 3, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 2, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 5, 3, 0, 1, 0, 0, 0, 0, 1, 3, 0, 3, 4, 0, 0, 2, 1, 0, 2, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 3, 5, 0, 0, 0, 0, 0, 1, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 3, 2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2, 0, 1, 0, 3, 1, 2, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 4, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 4, 1, 2, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 3, 0, 2, 2, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 2, 0, 0, 1, 3, 2, 2, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 4, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 2, 1, 1, 1, 0, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 2, 1, 1, 2, 3, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1};
        final double[] t = new double[arr.length];
        Arrays.fill(t, 1.0);

        final double mean = StatUtils.mean(arr);
        final double var = StatUtils.variance(arr);
        final double mom_alpha = (mean * mean) / (var - mean),
                     mom_beta = mean / (var - mean);

        final PointValuePair result = TLambdaPoissonSmoother.estimateAlphaBeta(arr, t);
//        System.err.format("MOM: %f\t%f\nML: %f\t%f\n", mom_alpha, mom_beta,
//                result.getPointRef()[0],
//                result.getPointRef()[1]);

        // Actual results are close to MOM result
        Assert.assertEquals(mom_alpha, result.getPointRef()[0], TOL);
        Assert.assertEquals(mom_beta, result.getPointRef()[1], TOL);

        final double[] tlSmoothed = TLambdaPoissonSmoother.smooth(arr, t);
        final double[] origSmoothed = EmpiricalBayesPoissonSmoother.smooth(arr);

        // this is a pretty loose tolerance...
        Assert.assertArrayEquals(origSmoothed, tlSmoothed, 0.1);
    }
}
