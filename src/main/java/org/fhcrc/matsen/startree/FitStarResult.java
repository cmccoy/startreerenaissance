package org.fhcrc.matsen.startree;

import java.util.List;
import java.util.Map;

/**
 * Serialized fit model
 * Created by cmccoy on 2/4/14.
 */
public class FitStarResult {
    private double[] branchLengths;
    private int degreesOfFreedom;
    private double logLikelihood, meanBranchLength;
    private Map<String, Double> independentParameters;
    private Map<String, Double> parameters;
    private List<Partition> partitions;

    public static class Partition {
        double[] Pmean;
        double[] Q;
        String modelName;
        double[] pi;
        Map<String, Double> parameters;
    }

    public double[] getBranchLengths() {
        return branchLengths;
    }

    public void setBranchLengths(double[] branchLengths) {
        this.branchLengths = branchLengths;
    }

    public int getDegreesOfFreedom() {
        return degreesOfFreedom;
    }

    public void setDegreesOfFreedom(int degreesOfFreedom) {
        this.degreesOfFreedom = degreesOfFreedom;
    }

    public double getLogLikelihood() {
        return logLikelihood;
    }

    public void setLogLikelihood(double logLikelihood) {
        this.logLikelihood = logLikelihood;
    }

    public double getMeanBranchLength() {
        return meanBranchLength;
    }

    public void setMeanBranchLength(double meanBranchLength) {
        this.meanBranchLength = meanBranchLength;
    }

    public Map<String, Double> getIndependentParameters() {
        return independentParameters;
    }

    public void setIndependentParameters(Map<String, Double> independentParameters) {
        this.independentParameters = independentParameters;
    }

    public Map<String, Double> getParameters() {
        return parameters;
    }

    public void setParameters(Map<String, Double> parameters) {
        this.parameters = parameters;
    }

    public List<Partition> getPartitions() {
        return partitions;
    }

    public void setPartitions(List<Partition> partitions) {
        this.partitions = partitions;
    }
}
