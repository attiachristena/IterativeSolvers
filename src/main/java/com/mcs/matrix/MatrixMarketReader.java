package com.mcs.matrix;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileReader;
import org.ejml.data.DMatrixSparseCSC;

public class MatrixMarketReader {

    public static DMatrixSparseCSC readMatrixMarketFile(String filePath) {
        DMatrixSparseCSC matrix = null;

        try (BufferedReader buffer = new BufferedReader(new FileReader(filePath))) {

            String line = buffer.readLine();

            // skip comment lines
            while (line != null && line.startsWith("%")) {
                line = buffer.readLine();
            }

            if (line == null) {
                throw new IOException("File vuoto o non valido");
            }

            String[] dimensions = line.trim().split("\\s+");
            int rows = Integer.parseInt(dimensions[0]);
            int cols = Integer.parseInt(dimensions[1]);
            int nonZeros = Integer.parseInt(dimensions[2]);

            matrix = new DMatrixSparseCSC(rows, cols, nonZeros);

            // lettura elementi
            while ((line = buffer.readLine()) != null) {
                if (line.startsWith("%") || line.trim().isEmpty()) continue;

                String[] values = line.trim().split("\\s+");

                int row = Integer.parseInt(values[0]) - 1;
                int col = Integer.parseInt(values[1]) - 1;
                double value = Double.parseDouble(values[2]);

                matrix.set(row, col, value);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        return matrix;
    }
}