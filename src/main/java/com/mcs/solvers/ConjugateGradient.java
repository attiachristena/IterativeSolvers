package com.mcs.solvers;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

import com.mcs.matrix_utils.MatrixValidator;
import com.mcs.result.SolverResult;

public class ConjugateGradient implements IterativeSolver{
    @Override
    public SolverResult solve(DMatrixSparseCSC A, DMatrixRMaj b, double tol, int maxIter) {

        // Initial checks
        if (!MatrixValidator.isSquare(A)) {
            throw new IllegalArgumentException("La matrice non è quadrata");
        }

        int n = A.numRows;
         
        DMatrixRMaj x = new DMatrixRMaj(n, 1);
        DMatrixRMaj r = new DMatrixRMaj(n, 1);
        DMatrixRMaj p = new DMatrixRMaj(n, 1);
        
        // x0 = null vector
        // r0 = b - A*x0 = b
        r.setTo(b);
        
        // p = r
        p.setTo(r);

        // Temporary vector: Ap = A*p
        DMatrixRMaj Ap = new DMatrixRMaj(n, 1);
        DMatrixRMaj r_new = new DMatrixRMaj(n, 1);

        // Initial relative norm
        double normB = NormOps_DDRM.normF(b);
        double normR = NormOps_DDRM.normF(r);
        double ratio = normR / normB;

        int iterations = 0;

        long startTime = System.currentTimeMillis();

        while (iterations < maxIter && ratio >= tol) {

            // Ap = A*p
            CommonOps_DSCC.mult(A, p, Ap);

            // rr = rᵀr
            double rr = CommonOps_DDRM.dot(r, r);

            // pAp = pᵀAp
            double pAp = CommonOps_DDRM.dot(p, Ap);

            // alpha = (rᵀr)/(pᵀAp)
            double alpha = rr / pAp;

            // x = x + alpha*p
            CommonOps_DDRM.addEquals(x, alpha, p);

            r_new.setTo(r);  // copio r in r_new

            // r_new = r - alpha*Ap
            CommonOps_DDRM.addEquals(r_new, -alpha, Ap);

            // rNrN = r_new * r_new
            double rNrN = CommonOps_DDRM.dot(r_new, r_new);

            // β = (r_new·r_new) / (r·r)
            double beta = rNrN / rr;

            // p = r_new + β*p
            CommonOps_DDRM.add(r_new, beta, p, p);

            r.setTo(r_new);

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
