package com.mcs.matrix_utils;
import org.ejml.data.DMatrixSparseCSC;

public class MatrixValidator {

    public static boolean isSymmetric(DMatrixSparseCSC A, double tol) {

        if (!isSquare(A)) {
            return false;
        }

        for (int col = 0; col < A.numCols; col++) {
            int start = A.col_idx[col];
            int end = A.col_idx[col + 1];

            for (int i = start; i < end; i++){
                int row = A.nz_rows[i];
                double value = A.nz_values[i]; // A[row, col]

                if (Math.abs(value - A.get(col, row)) > tol) 
                    return false;  
            }
        }

        return true;
    }

    public static boolean isSquare(DMatrixSparseCSC A) {
        return A.numRows == A.numCols;
    }

    public static boolean hasNonZeroDiagonal (DMatrixSparseCSC A){
        for (int i = 0; i < A.numRows; i++){
            if (A.get(i, i) == 0){
                return false;
            }
        }
        return true;
    }

    public static void validateMatrix(DMatrixSparseCSC A) throws IllegalArgumentException {
        if (A == null) {
            throw new IllegalArgumentException("La matrice non può essere null");
        }
        if (A.numRows <= 0 || A.numCols <= 0) {
            throw new IllegalArgumentException("La matrice deve avere dimensioni positive");
        }
    }

}
