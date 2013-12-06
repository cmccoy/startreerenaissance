package org.fhcrc.matsen.startree;

import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;
import com.google.common.base.Preconditions;

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

    public void print(final java.io.PrintStream ps) {
        DoubleMatrix2D[] matrices = new DoubleMatrix2D[] {
                conditionalNonsynonymous,
                conditionalSynonymous,
                unconditionalNonsynonymous,
                unconditionalSynonymous
        };
        String[] types = new String[] { "N", "S", "N", "S"};
        String[] conditions = new String[] { "C", "C", "U", "U"};
        ps.print("state");
        for(String condition : conditions) {
            for(String type : types) {
                for(int i = 0; i < conditionalNonsynonymous.columns(); i++) {
                    ps.print('\t');
                    ps.format("%s%s[%d]", condition, type, i + 1);
                }
            }
        }
        ps.print('\n');

        for(int row = 0; row < matrices[0].rows(); row++) {
            ps.format("%d", row + 1);
            for(final DoubleMatrix2D m : matrices) {
                for(int col = 0; col < matrices[0].columns(); col++) {
                    ps.print('\t');
                    ps.print(m.get(row, col));
                }
            }
            ps.print('\n');
        }
    }
}
