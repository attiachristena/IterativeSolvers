package com.mcs.solvers;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;

import com.mcs.result.SolverResult;

public interface IterativeSolver {
    
    public SolverResult solve(DMatrixSparseCSC A, DMatrixRMaj b, double tol, int maxIter);

}
