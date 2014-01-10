package org.fhcrc.matsen.startree;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceDistribution;
import org.apache.commons.math.linear.RealMatrix;

import java.util.List;
import java.util.Map;

/**
 * Created by cmccoy on 1/8/14.
 */
public class StarTreeRenaissanceResult {
    private StarTreeTraces unsmoothed, smoothed;
    private Map<String, List<TraceDistribution>> traceDistributions;

    // For serialization
    protected StarTreeRenaissanceResult() {

    }

    public StarTreeRenaissanceResult(final StarTreeTraces unsmoothed, final boolean sample, final int burnin) {
        this.unsmoothed = unsmoothed;
        this.smoothed = unsmoothed.getSmoothed(sample);
        this.traceDistributions = getTraceDistributions(smoothed, burnin);
    }

    public StarTreeTraces getUnsmoothed() {
        return unsmoothed;
    }

    public StarTreeTraces getSmoothed() {
        return smoothed;
    }

    public Map<String, List<TraceDistribution>> getTraceDistributions() {
        return traceDistributions;
    }

    private static Map<String, List<TraceDistribution>> getTraceDistributions(final StarTreeTraces starTraces, final int burnin) {
        final Map<String, List<TraceDistribution>> result = Maps.newHashMap();

        final RealMatrix dNdS = starTraces.getDNdSMatrix();
        final List<Trace> traces = Lists.newArrayListWithCapacity(dNdS.getColumnDimension() + 1);

        final Trace stateTrace = new Trace("state");
        for(final double d : starTraces.getState().toArray())
            stateTrace.add(d);
        traces.add(stateTrace);
        for(int col = 0; col < dNdS.getColumnDimension(); col++) {
            final double[] colValues = dNdS.getColumn(col);
            Trace t = new Trace(String.format(String.format("dNdS[%d]", col + 1)));
            for(final double d : colValues)
                t.add(d);
            Preconditions.checkState(t.getValuesSize() == colValues.length,
                    "Expected %d values, got %d.",
                    colValues.length,
                    t.getValuesSize());
            traces.add(t);
        }

        ArrayTraceList traceList = new ArrayTraceList("traces", traces, burnin);
        final List<TraceDistribution> traceDists = Lists.newArrayListWithCapacity(dNdS.getColumnDimension());
        for(int i = 1; i < traceList.getTraceCount(); i++) {
            Preconditions.checkState(traceList.getTraceName(i).startsWith("dNdS"));
            traceList.analyseTrace(i);
            traceDists.add(traceList.getDistributionStatistics(i));
        }

        result.put("dNdS", traceDists);

        return result;
    }
}
