package org.fhcrc.matsen.util;

import com.google.common.base.Preconditions;

import java.util.Arrays;

/**
 * Created by cmccoy on 1/6/14.
 */
public class Doubles {
    private Doubles() {
        throw new IllegalStateException("Do not instantiate.");
    }

    public static double[] constantArray(final int length, final double value) {
        Preconditions.checkArgument(length >= 0, "Positive length required (got: %s)", length);
        final double[] result = new double[length];
        Arrays.fill(result, value);
        return result;
    }

    public static double[] scalarMultiply(final double[] array, final double scalar) {
        Preconditions.checkNotNull(array);
        final double[] result = Arrays.copyOf(array, array.length);

        for(int i = 0; i < result.length; i++)
            result[i] *= scalar;

        return result;
    }

}
