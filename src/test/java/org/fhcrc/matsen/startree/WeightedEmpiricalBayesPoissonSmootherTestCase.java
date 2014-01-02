package org.fhcrc.matsen.startree;

import java.util.Arrays;

import dr.math.EmpiricalBayesPoissonSmoother;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by cmccoy on 12/16/13.
 */
@Deprecated
public class WeightedEmpiricalBayesPoissonSmootherTestCase {
    public static final double TOL = 1e-5;
    /**
     * Verify the unweighted case against the original Poisson smoothing code
     * @throws Exception
     */
    @Test
    public void testSmooth() throws Exception {
        final int n = 100;
        double[] arr = new double[n], w = new double[n];
        Arrays.fill(w, 1.0);
        Arrays.fill(arr, 0, n / 2, 20.0);
        Arrays.fill(arr, n / 2, n, 0.0);

        final double[] weightedSmoothed = WeightedEmpiricalBayesPoissonSmoother.smooth(arr, w);
        final double[] origSmoothed = EmpiricalBayesPoissonSmoother.smooth(arr);

        Assert.assertArrayEquals(origSmoothed, weightedSmoothed, TOL);
    }
    /**
     * Verify the unweighted case against the original Poisson smoothing code incorporating sampling
     * @throws Exception
     */
    @Test
    public void testSmoothSample() throws Exception {
        final int n = 100;
        double[] arr = new double[n], w = new double[n];
        Arrays.fill(w, 1.0);
        Arrays.fill(arr, 0, n / 2, 5.0);
        Arrays.fill(arr, n / 2, n, 1.0);

        final int seed = 1;
        dr.math.MathUtils.setSeed(seed);
        final double[] weightedSmoothed = WeightedEmpiricalBayesPoissonSmoother.smooth(arr, w, true);

        dr.math.MathUtils.setSeed(seed);
        final double[] origSmoothed = EmpiricalBayesPoissonSmoother.smoothWithSample(arr);

        Assert.assertArrayEquals(origSmoothed, weightedSmoothed, TOL);
    }
}
