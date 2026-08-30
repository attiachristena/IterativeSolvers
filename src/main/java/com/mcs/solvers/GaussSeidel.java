package com.mcs.solvers;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;

import com.mcs.matrix_utils.MatrixValidator;
import com.mcs.result.SolverResult;
import com.mcs.utils.UtilsOperations;

public class GaussSeidel implements IterativeSolver{

    @Override
    public SolverResult solve(DMatrixSparseCSC A, DMatrixRMaj b, double tol, int maxIter) {
        if (!MatrixValidator.isSquare(A)){
            throw new IllegalArgumentException("La matrice non è quadrata");
        }

        if (!MatrixValidator.hasNonZeroDiagonal(A)){
            throw new IllegalArgumentException("La matrice ha ha zeri sulla diagonale");
        }

        int n = A.numRows;

        DMatrixRMaj x = new DMatrixRMaj(n, 1);

        // Obtain references to internal arrays for more efficient access
        double[] xData = x.data;
        double[] bData = b.data;
        
        // Extract the diagonal
        double diagonal[] = UtilsOperations.extractDiagonal(A);

        long execTime = System.currentTimeMillis();

        int iterations = 0;
        double norm = 1.0;  

        // Iterate until convergence or max iterations
        while (iterations < maxIter && norm >= tol) {
            
            for (int i = 0; i < n; i++) {

                double sigma = 0.0; 

                int start = A.col_idx[i];
                int end = A.col_idx[i + 1];

                // Calculate the sum of products A[i][j] * x[j] for j != i
                for (int k = start; k < end; k++) {

                    int j = A.nz_rows[k];
                    double value = A.nz_values[k];

                    // For j != i, use the updated value of x[j] if j < i, otherwise use the previous value
                    if (j != i) {
                        sigma += value * xData[j]; 
                    }
                }

                // Update x[i] immediately using the Gauss-Seidel formula
                xData[i] = (bData[i] - sigma) / diagonal[i];
            }

            // Calculate residue
            norm = UtilsOperations.relativeResidue(A, x, b);
            iterations++;
        }

        double endTime = (System.currentTimeMillis() - execTime) / 1000.0;
        boolean converged = (norm < tol);

        return new SolverResult(x, iterations, norm, endTime, converged);
    
    }
    
}
