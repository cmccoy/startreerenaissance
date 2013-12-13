package org.fhcrc.matsen.startree;

import com.google.common.base.Preconditions;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.StatUtils;

import java.util.logging.Level;
import java.util.logging.Logger;


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
        for(int i = 0; i < m.getRowDimension(); i++) {
            result[i] = StatUtils.sum(m.getRow(i));
        }
        return result;
    }

}
