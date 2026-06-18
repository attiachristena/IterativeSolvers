package com.mcs.utils;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.sparse.csc.CommonOps_DSCC;

public class UtilsOperations {
    
    /*
    This class was created to accommodate any utility functions that might be needed in multiple parts 
    of the project, such as extracting the diagonal of a matrix or calculating the relative norm of the residue, 
    that are not directly tied to a particular solver.
    */

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
