package org.fhcrc.matsen.startree;

import dr.math.EmpiricalBayesPoissonSmoother;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.StatUtils;
import org.fhcrc.matsen.util.Doubles;
import org.junit.Assert;
import org.junit.Test;

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
public class TLambdaPoissonSmootherTestCase {
    public final double TOL = 3e-2;


    @Test
    public void testUnderdispersedSmooth() throws Exception {
        final double[] c = new double[]{ 10.0, 16.0, 25.0, 10.0, 13.0, 15.0, 18.0, 16.0, 18.0, 14.0, 19.0, 18.0, 17.0, 15.0, 19.0, 10.0, 17.0, 19.0, 13.0, 11.0, 10.0, 10.0, 14.0, 16.0, 7.0, 11.0, 15.0, 16.0, 15.0, 15.0, 14.0, 12.0, 19.0, 20.0, 12.0, 13.0, 21.0, 16.0, 15.0, 12.0, 9.0, 5.0, 12.0, 14.0, 14.0, 13.0, 15.0, 24.0, 12.0, 19.0, 15.0, 16.0, 16.0, 19.0, 11.0, 13.0, 17.0, 16.0, 12.0, 10.0, 16.0, 25.0, 14.0, 11.0, 20.0, 13.0, 16.0, 16.0, 14.0, 12.0, 23.0, 14.0, 15.0, 9.0, 11.0, 13.0, 19.0, 13.0, 11.0, 10.0, 15.0, 16.0, 11.0, 12.0, 15.0, 20.0, 14.0, 22.0, 11.0, 9.0, 13.0, 15.0, 24.0, 17.0, 19.0, 15.0, 17.0, 15.0, 15.0, 11.0 };
        final double[] t = Doubles.constantArray(c.length, 1.0);

        final double fitLambda = TLambdaPoissonSmoother.fitLambda(c, t);
        final double[] smoothed = TLambdaPoissonSmoother.smooth(c, t, false);
        final double[] smoothedBl = Doubles.scalarMultiply(TLambdaPoissonSmoother.smooth(c, Doubles.constantArray(c.length, 10.0), false), 10);

        for(int i = 0; i < smoothed.length; i++) {
            Assert.assertEquals(fitLambda, smoothed[i], 5e-2);
            Assert.assertEquals(fitLambda, smoothedBl[i], 5e-2);
        }
    }

    @Test
    public void testSmooth() throws Exception {
        final double[] arr = new double[] {0, 1, 3, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 4, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 2, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 1, 1, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 2, 0, 3, 1, 0, 0, 0, 0, 2, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 3, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 1, 1, 1, 0, 0, 2, 2, 3, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 6, 1, 0, 3, 1, 1, 3, 0, 2, 0, 1, 0, 0, 0, 2, 1, 4, 1, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 0, 0, 1, 1, 5, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 2, 1, 0, 2, 1, 0, 0, 0, 1, 0, 0, 2, 3, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 2, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 5, 3, 0, 1, 0, 0, 0, 0, 1, 3, 0, 3, 4, 0, 0, 2, 1, 0, 2, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 3, 5, 0, 0, 0, 0, 0, 1, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 3, 2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2, 0, 1, 0, 3, 1, 2, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 4, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 4, 1, 2, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 3, 0, 2, 2, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 2, 0, 0, 1, 3, 2, 2, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 4, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 2, 1, 1, 1, 0, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 2, 1, 1, 2, 3, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1};
        final double[] t = Doubles.constantArray(arr.length, 1.0);
        final double largeT = 10.0;

        final double mean = StatUtils.mean(arr);
        final double var = StatUtils.variance(arr);
        final double mom_alpha = (mean * mean) / (var - mean),
                     mom_beta = mean / (var - mean);

        final double[] result = TLambdaPoissonSmoother.estimateAlphaBeta(arr, t).getPointRef();
        final double[] resultBl = TLambdaPoissonSmoother.estimateAlphaBeta(arr, Doubles.constantArray(t.length, largeT)).getPointRef();

        // Larger branch length just changes Beta
        Assert.assertEquals(result[0], resultBl[0], 1e-2);
        Assert.assertEquals(result[1], resultBl[1] / largeT, 1e-2);

        // Actual results are close to MOM result
        Assert.assertEquals(mom_alpha, result[0], TOL);
        Assert.assertEquals(mom_beta, result[1], TOL);

        final double[] tlSmoothed = TLambdaPoissonSmoother.smooth(arr, t);
        final double[] tlSmoothedBl = TLambdaPoissonSmoother.smooth(arr, Doubles.constantArray(t.length, largeT));
        Assert.assertArrayEquals(tlSmoothed, Doubles.scalarMultiply(tlSmoothedBl, largeT), 1e2);
        final double[] origSmoothed = EmpiricalBayesPoissonSmoother.smooth(arr);

        // this is a pretty loose tolerance...
        Assert.assertArrayEquals(origSmoothed, tlSmoothed, 0.1);
    }


}
