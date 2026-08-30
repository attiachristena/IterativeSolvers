package com.mcs.utils;

import com.mcs.result.SolverResult;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Locale;

public class CSVExporter {

    private static final String OUTPUT_FOLDER = "results";
    private static final String OUTPUT_FILE = OUTPUT_FOLDER + "/results.csv";

    public static void createFile() throws IOException {

        File folder = new File(OUTPUT_FOLDER);
        if (!folder.exists()) {
            folder.mkdirs();
        }

        PrintWriter writer = new PrintWriter(new FileWriter(OUTPUT_FILE));

        writer.println("Matrix,Tolerance,Method,Iterations,ExecutionTime,RelativeError,RelativeResidual");

        writer.close();
    }

    public static void appendResult(String matrixName,
                                    double tolerance,
                                    String method,
                                    SolverResult result,
                                    double relativeError) throws IOException {

        PrintWriter writer = new PrintWriter(new FileWriter(OUTPUT_FILE, true));

        writer.printf(Locale.US, "%s,%e,%s,%d,%.6f,%e,%e%n",
        matrixName,
        tolerance,
        method,
        result.getIterations(),
        result.getExecutionTime(),
        relativeError,
        result.getRelativeResidual());

        writer.close();
    }
}