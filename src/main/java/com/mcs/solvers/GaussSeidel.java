package com.mcs.solvers;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;

import com.mcs.matrix.MatrixValidator;
import com.mcs.result.SolverResult;
import com.mcs.utils.UtilsOperations;

public class GaussSeidel implements IterativeSolver{

    // Funzionamento: simile al Jacobi, ma invece di calcolare tutti i nuovi valori di x in un nuovo vettore xNew, 
    // aggiorna direttamente il vettore x durante l'iterazione.
    // Questo significa che quando calcola x[i], utilizza i valori più aggiornati di x[j] per j < i, e i valori precedenti di x[j] per j > i. 
    // Questo può portare a una convergenza più rapida rispetto al Jacobi, soprattutto per matrici con forti interazioni tra le variabili.
    // Non viene fatta la suddivisione esplicita della matrice A in Lower Triangular e Upper Triangular e Diagonal, 
    // ma si sfrutta la struttura della matrice durante l'iterazione.

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

        // Ottengo i riferimenti agli array interni per un accesso più efficiente
        double[] xData = x.data;
        double[] bData = b.data;
        
        // Estraggo la diagonale
        double diagonal[] = UtilsOperations.extractDiagonal(A);

        long execTime = System.currentTimeMillis();

        int iterations = 0;
        double norm = 1.0;  

        // Iterazione fino a raggiungere la tolleranza o il numero massimo di iterazioni
        while (iterations < maxIter && norm >= tol) {
            
            for (int i = 0; i < n; i++) {

                double sigma = 0.0; 

                int start = A.col_idx[i];
                int end = A.col_idx[i + 1];

                // Calcolo la somma dei prodotti A[i][j] * x[j] per j != i
                for (int k = start; k < end; k++) {

                    int j = A.nz_rows[k];
                    double value = A.nz_values[k];

                    // Per j != i, uso il valore aggiornato di x[j] se j < i, altrimenti uso il valore precedente
                    if (j != i) {
                        sigma += value * xData[j]; 
                    }
                }

                // Aggiorno subito x[i] usando la formula di Gauss-Seidel
                xData[i] = (bData[i] - sigma) / diagonal[i];
            }

            // Calcolo residuo
            norm = UtilsOperations.relativeResidue(A, x, b);
            iterations++;
        }

        double endTime = (System.currentTimeMillis() - execTime) / 1000.0;
        return new SolverResult(x, iterations, norm, endTime, norm < tol);
    
    }
    
}
