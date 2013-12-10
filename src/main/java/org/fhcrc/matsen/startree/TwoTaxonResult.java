package org.fhcrc.matsen.startree;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.math.Functions;
import com.google.common.base.Preconditions;
import dr.math.EmpiricalBayesPoissonSmoother;


/**
 * Created with IntelliJ IDEA.
 * User: cmccoy
 * Date: 12/2/13
 * Time: 12:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class TwoTaxonResult implements java.io.Serializable {
    final DoubleMatrix2D conditionalNonsynonymous, conditionalSynonymous,
            unconditionalNonsynonymous, unconditionalSynonymous;

    public TwoTaxonResult(DoubleMatrix2D conditionalNonsynonymous, DoubleMatrix2D conditionalSynonymous, DoubleMatrix2D unconditionalNonsynonymous, DoubleMatrix2D unconditionalSynonymous) {
        Preconditions.checkArgument(conditionalNonsynonymous.rows() == conditionalSynonymous.rows());
        Preconditions.checkArgument(conditionalNonsynonymous.rows() == unconditionalNonsynonymous.rows());
        Preconditions.checkArgument(conditionalNonsynonymous.rows() == unconditionalSynonymous.rows());
        Preconditions.checkArgument(conditionalNonsynonymous.columns() == conditionalSynonymous.columns());
        Preconditions.checkArgument(conditionalNonsynonymous.columns() == unconditionalNonsynonymous.columns());
        Preconditions.checkArgument(conditionalNonsynonymous.columns() == unconditionalSynonymous.columns());

        this.conditionalNonsynonymous = conditionalNonsynonymous;
        this.conditionalSynonymous = conditionalSynonymous;
        this.unconditionalNonsynonymous = unconditionalNonsynonymous;
        this.unconditionalSynonymous = unconditionalSynonymous;
    }

    public TwoTaxonResult plus(final TwoTaxonResult other) {
        Preconditions.checkArgument(other.conditionalNonsynonymous.rows() == conditionalNonsynonymous.rows());
        Preconditions.checkArgument(other.conditionalNonsynonymous.columns() == conditionalNonsynonymous.columns());

        DoubleMatrix2D cn = conditionalNonsynonymous.copy(),
                cs = conditionalSynonymous.copy(),
                un = unconditionalNonsynonymous.copy(),
                us = unconditionalNonsynonymous.copy();
        cn.assign(other.conditionalNonsynonymous, Functions.plus);
        cs.assign(other.conditionalSynonymous, Functions.plus);
        un.assign(other.unconditionalNonsynonymous, Functions.plus);
        us.assign(other.unconditionalSynonymous, Functions.plus);

        return new TwoTaxonResult(cn, cs, un, us);
    }

    public DoubleMatrix2D getUnconditionalSynonymous() {
        return unconditionalSynonymous;
    }

    public DoubleMatrix2D getUnconditionalNonsynonymous() {
        return unconditionalNonsynonymous;
    }

    public DoubleMatrix2D getConditionalSynonymous() {
        return conditionalSynonymous;
    }

    public DoubleMatrix2D getConditionalNonsynonymous() { return conditionalNonsynonymous; }

    /**
     * Get a smoothed equivalent of this result. 
     *
     * This applies the Empirical Bayes smoothing of Lemey et. al. to each row of each matrix.
     */
    public TwoTaxonResult getSmoothed() {
        double[][] cn = conditionalNonsynonymous.toArray(),
                   cs = conditionalSynonymous.toArray(),
                   un = unconditionalNonsynonymous.toArray(),
                   us = unconditionalSynonymous.toArray();

        double[][][] arrays = new double[][][]{cn, cs, un, us};

        for(int i = 0; i < arrays.length; i++) {
            for(int j = 0; j < arrays[i].length; j++) {
                arrays[i][j] = EmpiricalBayesPoissonSmoother.smooth(arrays[i][j]);
            }
        }

        return new TwoTaxonResult(
            new DenseDoubleMatrix2D(cn),
            new DenseDoubleMatrix2D(cs),
            new DenseDoubleMatrix2D(un),
            new DenseDoubleMatrix2D(us));
    }

    /**
     * Apply the robust counting method of Lemey et. al. - this corresponds to Equation 1.
     */
    public DoubleMatrix2D getDNdSMatrix() {
        final DoubleMatrix2D result = conditionalNonsynonymous.like();

        for(int i = 0; i < conditionalNonsynonymous.rows(); i++) {
            for(int j = 0; j < conditionalNonsynonymous.columns(); j++) {
                final double d = (conditionalNonsynonymous.getQuick(i, j) * unconditionalSynonymous.getQuick(i, j)) /
                                 (conditionalSynonymous.getQuick(i, j) * unconditionalNonsynonymous.getQuick(i, j));
                result.setQuick(i, j, d);
            }
        }

        return result;
    }

    public void print(final java.io.PrintStream ps) {
        print(ps, true);
    }

    public void print(final java.io.PrintStream ps, final boolean dNdS) {
        DoubleMatrix2D[] matrices = new DoubleMatrix2D[] {
                conditionalNonsynonymous,
                conditionalSynonymous,
                unconditionalNonsynonymous,
                unconditionalSynonymous
        };
        final DoubleMatrix2D dndsMatrix = dNdS ? getDNdSMatrix() : null;
        String[] types = new String[] { "N", "S", "N", "S"};
        String[] conditions = new String[] { "C", "C", "U", "U"};
        ps.print("state");
        for(int m = 0; m < matrices.length; m++) {
            for(int i = 0; i < conditionalNonsynonymous.columns(); i++) {
                ps.print('\t');
                ps.format("%s%s[%d]", conditions[m], types[m], i + 1);
            }

        }
        if(dNdS) {
            for(int i = 0; i < conditionalNonsynonymous.columns(); i++)
                ps.format("\tdNdS[%d]", i + 1);
        }
        ps.print('\n');

        for(int row = 0; row < matrices[0].rows(); row++) {
            ps.format("%d", 10 * row);
            for(final DoubleMatrix2D m : matrices) {
                for(int col = 0; col < matrices[0].columns(); col++) {
                    ps.print('\t');
                    ps.print(m.get(row, col));
                }
            }
            if(dndsMatrix != null) {
                for(int col = 0; col < matrices[0].columns(); col++) {
                    ps.print('\t');
                    ps.print(dndsMatrix.getQuick(row, col));
                }
            }

            ps.print('\n');
        }
    }
}