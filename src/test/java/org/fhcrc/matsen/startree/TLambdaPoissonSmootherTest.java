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

//    @Test
//    public void testTempSmooth() throws Exception {
//        final double[] c = new double[]{ 10.0, 16.0, 25.0, 10.0, 13.0, 15.0, 18.0, 16.0, 18.0, 14.0, 19.0, 18.0, 17.0, 15.0, 19.0, 10.0, 17.0, 19.0, 13.0, 11.0, 10.0, 10.0, 14.0, 16.0, 7.0, 11.0, 15.0, 16.0, 15.0, 15.0, 14.0, 12.0, 19.0, 20.0, 12.0, 13.0, 21.0, 16.0, 15.0, 12.0, 9.0, 5.0, 12.0, 14.0, 14.0, 13.0, 15.0, 24.0, 12.0, 19.0, 15.0, 16.0, 16.0, 19.0, 11.0, 13.0, 17.0, 16.0, 12.0, 10.0, 16.0, 25.0, 14.0, 11.0, 20.0, 13.0, 16.0, 16.0, 14.0, 12.0, 23.0, 14.0, 15.0, 9.0, 11.0, 13.0, 19.0, 13.0, 11.0, 10.0, 15.0, 16.0, 11.0, 12.0, 15.0, 20.0, 14.0, 22.0, 11.0, 9.0, 13.0, 15.0, 24.0, 17.0, 19.0, 15.0, 17.0, 15.0, 15.0, 11.0 };
//        final double[] t = new double[c.length];
//        Arrays.fill(t, 9.918035);
//
//        try(java.io.FileWriter w = new java.io.FileWriter("test.csv")) {
//          w.write("alpha,beta,ll\n");
//          final double n = 50;
//          for(int i = 1; i <= n; i++) {
//              for(int j = 1; j <= n; j++) {
//                final double alpha =  i * 40.0 / n;
//                final double beta =  j * 40.0 / n;
//                final double ll = TLambdaPoissonSmoother.logLikelihood(alpha, beta, c, t);
//                w.write(String.format("%f,%f,%f\n", alpha, beta, ll));
//              }
//          }
//        }
//    }

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
