package org.fhcrc.matsen.startree;

import cern.colt.matrix.DoubleMatrix2D;

/**
 * Created with IntelliJ IDEA.
 * User: cmccoy
 * Date: 12/2/13
 * Time: 12:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class TwoTaxonResult {
    final DoubleMatrix2D conditionalNonsynonymous, conditionalSynonymous,
            unconditionalNonsynonymous, unconditionalSynonymous;

    public TwoTaxonResult(DoubleMatrix2D conditionalNonsynonymous, DoubleMatrix2D conditionalSynonymous, DoubleMatrix2D unconditionalNonsynonymous, DoubleMatrix2D unconditionalSynonymous) {
        this.conditionalNonsynonymous = conditionalNonsynonymous;
        this.conditionalSynonymous = conditionalSynonymous;
        this.unconditionalNonsynonymous = unconditionalNonsynonymous;
        this.unconditionalSynonymous = unconditionalSynonymous;
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

    public DoubleMatrix2D getConditionalNonsynonymous() {
        return conditionalNonsynonymous;
    }
}
