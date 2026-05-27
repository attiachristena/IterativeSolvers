package com.mcs.solvers;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;

import com.mcs.result.SolverResult;
import com.mcs.matrix.*;
import com.mcs.utils.*;


public class JacobiSolver implements IterativeSolver{

    @Override
    public SolverResult solve(DMatrixSparseCSC A, DMatrixRMaj b, double tol, int maxIter) {
                
        if (!MatrixValidator.isSquare(A)){
            throw new IllegalArgumentException("La matrice non è quadrata");
        }

        if (!MatrixValidator.hasNonZeroDiagonal(A)){
            throw new IllegalArgumentException("La matrice ha ha zeri sulla diagonale");
        }

        int n = A.numRows;

        // Inizializzo i vettori necessari
        DMatrixRMaj x = new DMatrixRMaj(n, 1);
        DMatrixRMaj xNew = new DMatrixRMaj(n, 1);

        // Ottengo i riferimenti agli array interni per un accesso più efficiente
        double[] xData = x.data;
        double[] xNewData = xNew.data;
        double[] bData = b.data;

        // Vettore per memorizzare la somma dei prodotti A[i][j] * x[j] per j != i
        double sum[] = new double[n];
        
        // Estraggo la diagonale
        double diagonal[] = UtilsOperations.extractDiagonal(A);

        long execTime = System.currentTimeMillis();

        int iterations = 0;
        double norm = 1.0;  

        while (iterations < maxIter && norm >= tol) {

            // Resetto il vettore sum per la nuova iterazione
            java.util.Arrays.fill(sum, 0.0); 
            
            // Calcolo la somma dei prodotti A[i][j] * x[j] per j != i per ogni riga i
            for (int col = 0; col < n; col++) {
                int start = A.col_idx[col];
                int end = A.col_idx[col + 1];

                // Per ogni elemento non zero nella colonna col, aggiorno la somma per la riga corrispondente
                for (int k = start; k < end; k++){
                    int row = A.nz_rows[k];
                    double value = A.nz_values[k]; 

                    if (row != col){
                        sum[row] += value * xData[col];
                    }
                }
            }

            // Aggiornamento Jacobi 
            for (int i = 0; i < n; i++){
                xNewData[i] = (bData[i] - sum[i]) / diagonal[i];
            }

            // Calcolo residuo
            norm = UtilsOperations.relativeResidue(A, xNew, b);
            x.setTo(xNew);
            iterations++;
        }

        double endTime = (System.currentTimeMillis() - execTime) / 1000.0;
        return new SolverResult(xNew, iterations, norm, endTime, norm < tol);
    }  
}
