package com.mcs.result;

import org.ejml.data.DMatrixRMaj;

public class SolverResult {

    private DMatrixRMaj solution;
    private int iterations;
    private double relativeResidual;
    private double executionTime;
    private boolean isConverged;

    public SolverResult(DMatrixRMaj solution, int iterations, double finalResidual, double executionTime, boolean isConverged) {
        this.solution = solution;
        this.iterations = iterations;
        this.relativeResidual = finalResidual;
        this.executionTime = executionTime;
        this.isConverged = isConverged;
    }

    public DMatrixRMaj getSolution(){
        return solution;
    }

    public double getExecutionTime() {
        return executionTime;
    }

    public int getIterations() {
        return iterations;
    }

    public double getRelativeResidual() {
        return relativeResidual;
    }

    public boolean isConverged() {
        return isConverged;
    }
    
}
