package com.mcs;
import com.mcs.matrix.MatrixMarketReader;
import com.mcs.result.SolverResult;
import com.mcs.solvers.ConjugateGradient;
import com.mcs.solvers.GaussSeidel;
import com.mcs.solvers.Gradient;
import com.mcs.solvers.JacobiSolver;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixRMaj;

public class Main {
    public static void main(String[] args){
        String filepath = "src/main/java/com/mcs/data/spa2.mtx";

        try{
        DMatrixSparseCSC matrix = MatrixMarketReader.readMatrixMarketFile(filepath);
        DMatrixRMaj b = new DMatrixRMaj(matrix.numRows, 1);
        // 1. xTrue = [1, 1, 1, ..., 1]
        DMatrixRMaj xTrue = new DMatrixRMaj(matrix.numRows, 1);
        for (int i = 0; i < matrix.numRows; i++) {
            xTrue.set(i, 0, 1.0);
        }

        // 2. b = A * xTrue
        org.ejml.sparse.csc.CommonOps_DSCC.mult(matrix, xTrue, b);

        if (matrix != null){
            System.out.println("Numero di colonne:" + matrix.numCols);
            System.out.println("Numero di righe:" + matrix.numRows);
            System.out.println("Numero di elementi non zero:" + matrix.nz_length);
        }

        JacobiSolver jacobi = new JacobiSolver();
        SolverResult result = jacobi.solve(matrix, b, 1e-6, 20000);

        GaussSeidel gaussSeidel = new GaussSeidel();
        SolverResult result2 = gaussSeidel.solve(matrix, b, 1e-6, 20000);

        Gradient gradient = new Gradient();
        SolverResult result3 = gradient.solve(matrix, b, 1e-6, 20000);

        ConjugateGradient conjugateGradient = new ConjugateGradient();
        SolverResult result4 = conjugateGradient.solve(matrix, b, 1e-6, 20000);

        System.out.println("Risultati per Jacobi");
        System.out.println("Converso: " + result.isConverged());
        System.out.println("Iterazioni: " + result.getIterations());
        System.out.println("Residuo: " + result.getRelativeError());
        System.out.println("Tempo: " + result.getExecutionTime() + "s");

        System.out.println("Risultati per Gauss-Seidel");
        System.out.println("Converso: " + result2.isConverged());
        System.out.println("Iterazioni: " + result2.getIterations());
        System.out.println("Residuo: " + result2.getRelativeError());
        System.out.println("Tempo: " + result2.getExecutionTime() + "s");

        System.out.println("Risultati per Gradiente");
        System.out.println("Converso: " + result3.isConverged());
        System.out.println("Iterazioni: " + result3.getIterations());
        System.out.println("Residuo: " + result3.getRelativeError());
        System.out.println("Tempo: " + result3.getExecutionTime() + "s");

        System.out.println("Risultati per Gradiente Coniugato");
        System.out.println("Converso: " + result4.isConverged());
        System.out.println("Iterazioni: " + result4.getIterations());
        System.out.println("Residuo: " + result4.getRelativeError());
        System.out.println("Tempo: " + result4.getExecutionTime() + "s");

    } catch (Exception e) {
        System.out.println("Si è verificato un errore durante la lettura del file:");
        e.printStackTrace();
    }
}
    
}
