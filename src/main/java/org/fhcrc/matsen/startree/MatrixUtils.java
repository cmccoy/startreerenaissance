package org.fhcrc.matsen.startree;

import com.google.common.base.Preconditions;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.BlockRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.StatUtils;


/**
 * Created with IntelliJ IDEA.
 * User: cmccoy
 * Date: 12/2/13
 * Time: 12:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class MatrixUtils {

    private MatrixUtils() {
        Preconditions.checkState(false);
    }

    public static double[] rowSums(final RealMatrix m) {
        double[] result = new double[m.getRowDimension()];
        for (int i = 0; i < m.getRowDimension(); i++) {
            result[i] = StatUtils.sum(m.getRow(i));
        }
        return result;
    }


    public static RealMatrix createRealMatrix(final int rows, final int cols, final double value) {
        RealMatrix result = new BlockRealMatrix(rows, cols);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.setEntry(i, j, value);
            }
        }

        return result;
    }

    public static RealVector createRealVector(final int length, final double value) {
        RealVector result = new ArrayRealVector(length);
        result.set(value);
        return result;
    }
}
