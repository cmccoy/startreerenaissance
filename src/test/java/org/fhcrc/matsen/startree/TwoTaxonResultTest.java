package org.fhcrc.matsen.startree;

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.BlockRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.junit.Test;
import org.junit.Assert;

/**
 * Created by cmccoy on 12/16/13.
 */
public class TwoTaxonResultTest {
    public static final double TOL = 1e-6;
    private static RealMatrix createRealMatrix(final int rows, final int cols, final double value) {
        RealMatrix result = new BlockRealMatrix(rows, cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.setEntry(i, j, value);
            }
        }

        return result;
    }

    private static RealVector createRealVector(final int length, final double value) {
        RealVector result = new ArrayRealVector(length);
        result.set(value);
        return result;
    }

    private void checkSum(final RealMatrix m1, final RealMatrix m2, final RealMatrix result) {
        Assert.assertEquals(m1.getColumnDimension(), m2.getColumnDimension());
        Assert.assertEquals(m1.getRowDimension(), m2.getRowDimension());
        Assert.assertEquals(m1.getColumnDimension(), result.getColumnDimension());
        Assert.assertEquals(m1.getRowDimension(), result.getRowDimension());
        for(int i = 0; i < m1.getRowDimension(); i++) {
            for(int j = 0; j < m1.getColumnDimension(); j++) {
                Assert.assertEquals(m1.getEntry(i, j) + m2.getEntry(i, j), result.getEntry(i, j), TOL);
            }
        }
    }

    @Test
    public void testPlus() throws Exception {
        final int r = 4, c = 12;
        final RealMatrix cn1 = createRealMatrix(r, c, 1.0),
                cn2 = createRealMatrix(r, c, 14.2),

                un1 = createRealMatrix(r, c, 12.1),
                un2 = createRealMatrix(r, c, 5.42),

                cs1 = createRealMatrix(r, c, 0.7),
                cs2 = createRealMatrix(r, c, 2.74),

                us1 = createRealMatrix(r, c, 7.7),
                us2 = createRealMatrix(r, c, 3.14);

        final double[] stateArr = new double[] {1, 20, 40, 80};
        final RealVector state = new ArrayRealVector(stateArr);

        final RealVector cov1 = createRealVector(c, 1.0);
        final RealVector cov2 = createRealVector(c, 1.0);

        final TwoTaxonResult t1 = new TwoTaxonResult(state, cn1, cs1, un1, us1, cov1.toArray()),
                             t2 = new TwoTaxonResult(state, cn2, cs2, un2, us2, cov2.toArray());

        final TwoTaxonResult summed = t1.plus(t2);

        checkSum(cn1, cn2, summed.getConditionalNonsynonymous());
        checkSum(un1, un2, summed.getUnconditionalNonsynonymous());
        checkSum(cs1, cs2, summed.getConditionalSynonymous());
        checkSum(us1, us2, summed.getUnconditionalSynonymous());

        final RealVector totalCov = cov1.add(cov2);
        Assert.assertArrayEquals(totalCov.getData(), summed.getCoverage().getData(), TOL);
    }
}