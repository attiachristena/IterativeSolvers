package com.mcs.solvers;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;

import com.mcs.matrix_utils.*;
import com.mcs.result.SolverResult;
import com.mcs.utils.*;


public class JacobiSolver implements IterativeSolver{

    @Override
    public SolverResult solve(DMatrixSparseCSC A, DMatrixRMaj b, double tol, int maxIter) {
                
        if (!MatrixValidator.isSquare(A)){
            throw new IllegalArgumentException("La matrice non è quadrata");
        }

        if (!MatrixValidator.hasNonZeroDiagonal(A)){
            throw new IllegalArgumentException("La matrice ha zeri sulla diagonale");
        }

        int n = A.numRows;

        // Inizialize solution vectors
        DMatrixRMaj x = new DMatrixRMaj(n, 1);
        DMatrixRMaj xNew = new DMatrixRMaj(n, 1);

        // Obtain references to internal arrays for more efficient access
        double[] xData = x.data;
        double[] xNewData = xNew.data;
        double[] bData = b.data;

        // Vector for storing the sum of products A[i][j] * x[j] for j != i
        double sum[] = new double[n];
        
        // Extract the diagonal
        double diagonal[] = UtilsOperations.extractDiagonal(A);

        long execTime = System.currentTimeMillis();

        int iterations = 0;
        double norm = 1.0;  

        while (iterations < maxIter && norm >= tol) {

            // Reset the sum vector for the new iteration
            java.util.Arrays.fill(sum, 0.0); 
            
            // Calculate the sum of products A[i][j] * x[j] for j != i for each row i
            for (int col = 0; col < n; col++) {
                int start = A.col_idx[col];
                int end = A.col_idx[col + 1];

                // For each non-zero element in column col, update the sum for the corresponding row
                for (int k = start; k < end; k++){
                    int row = A.nz_rows[k];
                    double value = A.nz_values[k]; 

                    if (row != col){
                        sum[row] += value * xData[col];
                    }
                }
            }

            // Update Jacobi
            for (int i = 0; i < n; i++){
                xNewData[i] = (bData[i] - sum[i]) / diagonal[i];
            }

            // Calculate residue
            norm = UtilsOperations.relativeResidue(A, xNew, b);
            x.setTo(xNew);
            iterations++;
        }

        double endTime = (System.currentTimeMillis() - execTime) / 1000.0;
        boolean converged = (norm < tol);
        return new SolverResult(xNew, iterations, norm, endTime, converged);
    }  
}
