package com.mcs.matrix_utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileReader;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.ops.DConvertMatrixStruct;

public class MatrixMarketReader {

    public static DMatrixSparseCSC readMatrixMarketFile(String filePath) {
        DMatrixSparseCSC matrix = null;

        try (BufferedReader buffer = new BufferedReader(new FileReader(filePath))) {

            String line = buffer.readLine();

            // Skip comment lines starting with '%'
            while (line != null && line.startsWith("%")) {
                line = buffer.readLine();
            }

            if (line == null) {
                throw new IOException("File vuoto o non valido");
            }

            // Read matrix dimensions and number of non-zero entries
            String[] dimensions = line.trim().split("\\s+"); 
            int rows = Integer.parseInt(dimensions[0]);
            int cols = Integer.parseInt(dimensions[1]);
            int nonZeros = Integer.parseInt(dimensions[2]);
            
            // Creates a triplet format matrix, which is more efficient for initial construction of the sparse matrix.
            DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(rows, cols, nonZeros); 
            
            
            while ((line = buffer.readLine()) != null) {
                if (line.startsWith("%") || line.trim().isEmpty()) continue;

                String[] values = line.trim().split("\\s+");

                int row = Integer.parseInt(values[0]) - 1;
                int col = Integer.parseInt(values[1]) - 1;
                double value = Double.parseDouble(values[2]);

                triplet.addItem(row, col, value);
            }

            // Converts the triplet format to a compressed sparse column (CSC) format, 
            // which is more efficient for matrix operations.
            matrix = new DMatrixSparseCSC(rows, cols, nonZeros);
            DConvertMatrixStruct.convert(triplet, matrix);

        } catch (IOException e) {
            e.printStackTrace();
        }

        return matrix;
    }
}