package org.fhcrc.matsen.startree;

import dr.math.EmpiricalBayesPoissonSmoother;
import org.apache.commons.math3.optim.PointValuePair;
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
    public final double TOL = 1e-5;

    @Test
    public void testSmooth() throws Exception {
        final double[] arr = new double[] {0, 1, 3, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 4, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 2, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 1, 1, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 2, 0, 3, 1, 0, 0, 0, 0, 2, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 3, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 1, 1, 1, 0, 0, 2, 2, 3, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 6, 1, 0, 3, 1, 1, 3, 0, 2, 0, 1, 0, 0, 0, 2, 1, 4, 1, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 0, 0, 1, 1, 5, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 2, 1, 0, 2, 1, 0, 0, 0, 1, 0, 0, 2, 3, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 2, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 5, 3, 0, 1, 0, 0, 0, 0, 1, 3, 0, 3, 4, 0, 0, 2, 1, 0, 2, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 3, 5, 0, 0, 0, 0, 0, 1, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 3, 2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2, 0, 1, 0, 3, 1, 2, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 4, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 3, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 4, 1, 2, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 3, 0, 2, 2, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 2, 0, 0, 1, 3, 2, 2, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 4, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 2, 1, 1, 1, 0, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 2, 1, 1, 2, 3, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1};
        final double[] t = new double[arr.length];
        Arrays.fill(t, 1.0);

        final PointValuePair result = TLambdaPoissonSmoother.estimateAlphaBeta(arr, t);
        // Actual results are just close
        Assert.assertEquals(1.0, result.getPointRef()[0], 0.2);
        Assert.assertEquals(2.0, result.getPointRef()[1], 0.2);

        final double[] tlSmoothed = TLambdaPoissonSmoother.smooth(arr, t);
        final double[] origSmoothed = EmpiricalBayesPoissonSmoother.smooth(arr);

        // this is a pretty loose tolerance...
        Assert.assertArrayEquals(origSmoothed, tlSmoothed, 0.1);
    }
}
