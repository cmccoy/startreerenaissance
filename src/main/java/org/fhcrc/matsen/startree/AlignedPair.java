package org.fhcrc.matsen.startree;

/**
 * Created by cmccoy on 12/5/13.
 */
public class AlignedPair {
    private final int queryPosition;
    private final int referencePosition;

    public AlignedPair(int queryPosition, int referencePosition) {
        this.queryPosition = queryPosition;
        this.referencePosition = referencePosition;
    }

    public int getQueryPosition() {
        return queryPosition;
    }

    public int getReferencePosition() {
        return referencePosition;
    }

    public boolean consumesQuery() {
        return queryPosition >= 0;
    }

    public boolean consumesReference() {
        return referencePosition >= 0;
    }
}
