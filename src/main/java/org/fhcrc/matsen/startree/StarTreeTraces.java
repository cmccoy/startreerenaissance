package org.fhcrc.matsen.startree;

import com.google.common.base.Preconditions;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.BlockRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.StatUtils;

import java.util.logging.Logger;


/**
 * Created with IntelliJ IDEA.
 * User: cmccoy
 * Date: 12/2/13
 * Time: 12:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class StarTreeTraces implements java.io.Serializable {
    private static final long serialVersionUID = 1;

    private static final Logger logger = Logger.getLogger("org.fhcrc.matsen.startree");

    final RealMatrix conditionalNonsynonymous, conditionalSynonymous,
            unconditionalNonsynonymous, unconditionalSynonymous, dNdS, totalBranchLength;
    final RealVector state, coverage;

    private static void checkDimensions(final RealMatrix a, final RealMatrix b) {
        Preconditions.checkArgument(a.getRowDimension() == b.getRowDimension());
        Preconditions.checkArgument(a.getColumnDimension() == b.getColumnDimension());
    }

    public StarTreeTraces(RealVector state, RealMatrix conditionalNonsynonymous, RealMatrix conditionalSynonymous, RealMatrix unconditionalNonsynonymous, RealMatrix unconditionalSynonymous, double[] coverage, RealMatrix totalBranchLength) {
        checkDimensions(conditionalNonsynonymous, conditionalSynonymous);
        checkDimensions(conditionalNonsynonymous, unconditionalNonsynonymous);
        checkDimensions(conditionalNonsynonymous, unconditionalSynonymous);
        Preconditions.checkArgument(conditionalNonsynonymous.getRowDimension() == state.getDimension());
        Preconditions.checkArgument(conditionalNonsynonymous.getColumnDimension() == coverage.length);
        checkDimensions(conditionalNonsynonymous, totalBranchLength);

        this.conditionalNonsynonymous = conditionalNonsynonymous;
        this.conditionalSynonymous = conditionalSynonymous;
        this.unconditionalNonsynonymous = unconditionalNonsynonymous;
        this.unconditionalSynonymous = unconditionalSynonymous;
        this.state = state;
        this.coverage = new ArrayRealVector(coverage);
        this.totalBranchLength = totalBranchLength;
        this.dNdS = computeDNdSMatrix();
    }

    public StarTreeTraces plus(final StarTreeTraces other) {
        Preconditions.checkArgument(other.conditionalNonsynonymous.getRowDimension() == conditionalNonsynonymous.getRowDimension());
        Preconditions.checkArgument(other.conditionalNonsynonymous.getColumnDimension() == conditionalNonsynonymous.getColumnDimension());

        RealMatrix cn = conditionalNonsynonymous.add(other.conditionalNonsynonymous),
                cs = conditionalSynonymous.add(other.conditionalSynonymous),
                un = unconditionalNonsynonymous.add(other.unconditionalNonsynonymous),
                us = unconditionalSynonymous.add(other.unconditionalSynonymous);

        return new StarTreeTraces(state, cn, cs, un, us, coverage.add(other.coverage).toArray(),
                totalBranchLength.add(other.totalBranchLength));
    }

    public RealMatrix getUnconditionalSynonymous() {
        return unconditionalSynonymous;
    }

    public RealMatrix getUnconditionalNonsynonymous() {
        return unconditionalNonsynonymous;
    }

    public RealMatrix getConditionalSynonymous() {
        return conditionalSynonymous;
    }

    public RealMatrix getConditionalNonsynonymous() {
        return conditionalNonsynonymous;
    }

    public RealVector getCoverage() {
        return coverage;
    }

    public RealMatrix getTotalBranchLength() {
        return totalBranchLength;
    }

    /**
     * Get a smoothed equivalent of this result.
     * <p/>
     * This applies the Empirical Bayes smoothing of Lemey et. al. to each row of each matrix.
     */
    public StarTreeTraces getSmoothed(boolean sample) {
        double[][] cn = conditionalNonsynonymous.getData(),
                cs = conditionalSynonymous.getData(),
                un = unconditionalNonsynonymous.getData(),
                us = unconditionalSynonymous.getData(),
                bl = totalBranchLength.getData();

        double[][][] arrays = new double[][][]{cn, cs, un, us};

        for (int arr = 0; arr < arrays.length; arr++) {
            for (int row = 0; row < arrays[arr].length; row++) {
                if (StatUtils.sum(arrays[arr][row]) == 0.0) {
                    logger.warning(String.format("No counts observed at array %d row %d", arr, row));
                }
                arrays[arr][row] = TLambdaPoissonSmoother.smooth(arrays[arr][row], bl[row], sample);
            }
        }

        return new StarTreeTraces(
                state.copy(),
                new BlockRealMatrix(cn),
                new BlockRealMatrix(cs),
                new BlockRealMatrix(un),
                new BlockRealMatrix(us),
                coverage.toArray(),
                totalBranchLength.copy());
    }

    public RealMatrix getDNdSMatrix() {
        return dNdS;
    }

    /**
     * Apply the robust counting method of Lemey et. al. - this corresponds to Equation 1.
     */
    private RealMatrix computeDNdSMatrix() {
        final RealMatrix result = conditionalNonsynonymous.createMatrix(conditionalNonsynonymous.getRowDimension(),
                conditionalNonsynonymous.getColumnDimension());

        for (int i = 0; i < conditionalNonsynonymous.getRowDimension(); i++) {
            for (int j = 0; j < conditionalNonsynonymous.getColumnDimension(); j++) {
                final double cn = conditionalNonsynonymous.getEntry(i, j),
                        un = unconditionalNonsynonymous.getEntry(i, j),
                        cs = conditionalSynonymous.getEntry(i, j),
                        us = unconditionalSynonymous.getEntry(i, j);
                final double d = (cn / un) / (cs / us);
                if (Double.isInfinite(d)) {
                    //logger.warning(String.format("Infinite dNdS for (%d, %d): cn=%f un=%f cs=%f us=%f - using 1.0", i, j, cn, un, cs, us));
                    result.setEntry(i, j, 1.0);
                } else {
                    result.setEntry(i, j, d);
                }
            }
        }

        return result;
    }


    public void print(final java.io.PrintStream ps) {
        print(ps, true);
    }

    public void print(final java.io.PrintStream ps, final boolean dNdS) {
        final com.google.common.base.Joiner joiner = com.google.common.base.Joiner.on('\t');
        RealMatrix[] matrices = new RealMatrix[]{
                conditionalNonsynonymous,
                conditionalSynonymous,
                unconditionalNonsynonymous,
                unconditionalSynonymous
        };

        ps.format("# Coverage: %s\n", java.util.Arrays.toString(coverage.toArray()));

        final RealMatrix dndsMatrix = dNdS ? getDNdSMatrix() : null;
        String[] types = new String[]{"N", "S", "N", "S"};
        String[] conditions = new String[]{"C", "C", "U", "U"};
        ps.print("state");
        for (int m = 0; m < matrices.length; m++) {
            for (int i = 0; i < conditionalNonsynonymous.getColumnDimension(); i++) {
                ps.print('\t');
                ps.format("%s%s[%d]", conditions[m], types[m], i + 1);
            }

        }

        ps.print(joiner.join("", "c_N", "c_S", "u_N", "u_S", "dNdS"));
        if (dNdS) {
            for (int i = 0; i < conditionalNonsynonymous.getColumnDimension(); i++)
                ps.format("\tdNdS[%d]", i + 1);
        }
        ps.print('\n');

        for (int row = 0; row < matrices[0].getRowDimension(); row++) {
            // Same order as matrices
            double[] sums = new double[]{0, 0, 0, 0};
            ps.print(state.getEntry(row));
            for (int i = 0; i < matrices.length; i++) {
                for (int col = 0; col < matrices[0].getColumnDimension(); col++) {
                    ps.print('\t');
                    ps.print(matrices[i].getEntry(row, col));
                    sums[i] += matrices[i].getEntry(row, col);
                }
            }

            // cn, cs, un, us, dnds
            ps.print(joiner.join("", sums[0], sums[1], sums[2], sums[3],
                    (sums[0] / sums[2]) / (sums[1] / sums[3])));

            if (dndsMatrix != null) {
                for (int col = 0; col < matrices[0].getColumnDimension(); col++) {
                    ps.print('\t');
                    ps.print(dndsMatrix.getEntry(row, col));
                }
            }

            ps.print('\n');
        }
    }

    public RealVector getState() {
        return state;
    }
}
