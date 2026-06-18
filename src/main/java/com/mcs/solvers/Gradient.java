package com.mcs.solvers;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.mcs.matrix.MatrixValidator;
import com.mcs.result.SolverResult;

public class Gradient implements IterativeSolver {

    @Override
    public SolverResult solve(DMatrixSparseCSC A, DMatrixRMaj b, double tol, int maxIter) {

        // Initial checks
        if (!MatrixValidator.isSquare(A)) {
            throw new IllegalArgumentException("La matrice non è quadrata");
        }

        int n = A.numRows;
         
        DMatrixRMaj x = new DMatrixRMaj(n, 1);
        DMatrixRMaj r = new DMatrixRMaj(n, 1);
        
        // x0 = null vector
        // r0 = b - A*x0 = b
        r.setTo(b);
        
        // Temporary vector: Ar = A*r
        DMatrixRMaj Ar = new DMatrixRMaj(n, 1);

        // Initial relative norm
        double normB = NormOps_DDRM.normF(b);
        double normR = NormOps_DDRM.normF(r);
        double ratio = normR / normB;

        int iterations = 0;

        long startTime = System.currentTimeMillis();

        while (iterations < maxIter && ratio >= tol) {

            // Ar = A*r
            CommonOps_DSCC.mult(A, r, Ar);

            // rr = rᵀr
            double rr = CommonOps_DDRM.dot(r, r);

            // rAr = rᵀAr
            double rAr = CommonOps_DDRM.dot(r, Ar);

            // alpha = (rᵀr)/(rᵀAr)
            double alpha = rr / rAr;

            // x = x + alpha*r
            CommonOps_DDRM.addEquals(x, alpha, r);

            // r = r - alpha*Ar
            CommonOps_DDRM.addEquals(r, -alpha, Ar);

            // Update relative residual
            normR = NormOps_DDRM.normF(r);
            ratio = normR / normB;

            iterations++;
        }

        double execTime = (System.currentTimeMillis() - startTime) / 1000.0;
        boolean converged = (ratio < tol);

        return new SolverResult(x, iterations, ratio, execTime, converged);
    }
}