package com.mcs.utils;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

public class UtilsOperations {
    // Questa classe è stata creata per ospitare eventuali funzioni di utilità che 
    // potrebbero essere necessarie in più parti del progetto, come ad esempio funzioni per il calcolo dell'errore relativo, 
    // la normalizzazione dei vettori, o altre operazioni comuni che non sono specifiche di un particolare solver.

    // ==== MATRIX OPERATIONS ====
    public static double[] extractDiagonal(DMatrixSparseCSC A){
        int n = A.numRows;
        double[] diagonal = new double[n];

        for (int col = 0; col < n; col++){
            int start = A.col_idx[col];
            int end = A.col_idx[col + 1];

            for (int k = start; k < end; k++){
                int row = A.nz_rows[k];
                double value = A.nz_values[k];

                if (row == col){
                    diagonal[col] = value;
                    break;
                }
            }
        }

        return diagonal;
    }

    // Check se la matrice è diagonalmente dominante
    // public static boolean isDiagonallyDominant(DMatrixSparseCSC A){

    // === SOLVER OPERATIONS ===
    public static double relativeResidue (DMatrixSparseCSC A, DMatrixRMaj x, DMatrixRMaj b){
        
        DMatrixRMaj residue = new DMatrixRMaj(A.numRows, 1);
        CommonOps_DSCC.mult(A, x, residue);
        CommonOps_DDRM.subtractEquals(residue, b);

        double norm = NormOps_DDRM.normF(residue) / NormOps_DDRM.normF(b);
        return norm;

    }

    public static double relativeError (DMatrixRMaj xExact, DMatrixRMaj xComputed){
        DMatrixRMaj diff = new DMatrixRMaj(xExact.numRows, 1);
        diff.setTo(xExact);

        CommonOps_DDRM.subtractEquals(diff, xComputed);
        return NormOps_DDRM.normF(diff) / NormOps_DDRM.normF(xExact);
    }
}
