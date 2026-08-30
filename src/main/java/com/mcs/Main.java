package com.mcs;

import com.mcs.matrix_utils.MatrixMarketReader;
import com.mcs.result.SolverResult;
import com.mcs.solvers.*;
import com.mcs.utils.UtilsOperations;
import com.mcs.utils.CSVExporter;


import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixRMaj;
import org.ejml.sparse.csc.CommonOps_DSCC;

import java.util.Scanner;

public class Main {

    public static void main(String[] args) {

        Scanner sc = new Scanner(System.in);

        System.out.println("==================================");
        System.out.println(" SOLVER SISTEMI LINEARI ");
        System.out.println("==================================");
        System.out.println("Scegli modalità:");
        System.out.println("1 - Esecuzione automatica (batch)");
        System.out.println("2 - Modalità interattiva");
        System.out.print("Scelta: ");

        int mode = sc.nextInt();

        if (mode == 1) {
            runBatch();
        } else if (mode == 2) {
            runInteractive(sc);
        } else {
            System.out.println("Scelta non valida");
        }
    }

    // ********************************************************
    // BATCH MODE 
    // ********************************************************
    private static void runBatch() {

        try {
        CSVExporter.createFile();
        } catch (Exception e) {
            System.out.println("Errore nella creazione del file CSV.");
            e.printStackTrace();
            return;
        }

        String[] matrices = {
            "spa1",
            "spa2",
            "vem1",
            "vem2"
        };

        double[] tolerances = {1e-4, 1e-6, 1e-8, 1e-10};

        IterativeSolver[] solvers = {
                new JacobiSolver(),
                new GaussSeidel(),
                new Gradient(),
                new ConjugateGradient()
        };

        String[] names = {
                "Jacobi",
                "Gauss-Seidel",
                "Gradient",
                "Conjugate Gradient"
        };

        for (String file : matrices) {

            System.out.println("\n=====================================");
            System.out.println("MATRICE: " + file);
            System.out.println("=====================================");

            DMatrixSparseCSC A =
            MatrixMarketReader.readMatrixMarketFile(
                "src/main/java/com/mcs/data/" + file + ".mtx"
            );
            
            int n = A.numRows;

            // xTrue = 1
            DMatrixRMaj xTrue = new DMatrixRMaj(n, 1);
            for (int i = 0; i < n; i++) {
                xTrue.set(i, 0, 1.0);
            }

            // b = A xTrue
            DMatrixRMaj b = new DMatrixRMaj(n, 1);
            CommonOps_DSCC.mult(A, xTrue, b);

            System.out.printf("%-18s %-10s %-10s %-15s %-15s %-10s\n",
                    "Metodo", "Tol", "Iter", "Residuo", "Errore", "Time");

            for (double tol : tolerances) {

                for (int i = 0; i < solvers.length; i++) {

                    SolverResult result =
                            solvers[i].solve(A, b, tol, 20000);

                    double error =
                            UtilsOperations.relativeError(
                                xTrue,
                                result.getSolution()
                            );

                            try {
                                CSVExporter.appendResult(
                                        file.replace("src/main/java/com/mcs/data/", "")
                                            .replace(".mtx", ""),
                                        tol,
                                        names[i],
                                        result,
                                        error
                                );
                            } catch (Exception e) {
                                e.printStackTrace();
}

                    System.out.printf(
                            "%-18s %-10.0e %-10d %-15.3e %-15.3e %-10.4f\n",
                            names[i],
                            tol,
                            result.getIterations(),
                            result.getRelativeResidual(),
                            error,
                            result.getExecutionTime()
                    );
                }

                System.out.println("------------------------------------------------------");
            }
        }
    }

    // **************************
    // INTERACTIVE MODE 
    // **************************
    private static void runInteractive(java.util.Scanner sc) {

    System.out.println("\n==================================");
    System.out.println("      MODALITÀ INTERATTIVA");
    System.out.println("==================================");

    // **************************
    // SELEZIONE MATRICE
    // **************************
    String[] matrices = {
            "spa1.mtx",
            "spa2.mtx",
            "vem1.mtx",
            "vem2.mtx"
    };

    System.out.println("\nSeleziona matrice:");
    for (int i = 0; i < matrices.length; i++) {
        System.out.println((i + 1) + " - " + matrices[i]);
    }

    System.out.print("Scelta: ");
    int m = sc.nextInt();
    String path = "src/main/java/com/mcs/data/" + matrices[m - 1];

    // **************************
    // SELEZIONE TOLLERANZA
    // **************************
    double[] tolerances = {1e-4, 1e-6, 1e-8, 1e-10};

    System.out.println("\nSeleziona tolleranza:");
    for (int i = 0; i < tolerances.length; i++) {
        System.out.println((i + 1) + " - " + tolerances[i]);
    }

    System.out.print("Scelta: ");
    int t = sc.nextInt();
    double tol = tolerances[t - 1];

    // **************************
    // MAX ITER
    // **************************
    System.out.print("\nInserisci max iterazioni (default 20000): ");
    int maxIter = sc.nextInt();

    // **************************
    // SELEZIONE DEL METODO
    // **************************
    System.out.println("\nSeleziona metodo:");
    System.out.println("1 - Jacobi");
    System.out.println("2 - Gauss-Seidel");
    System.out.println("3 - Gradiente");
    System.out.println("4 - Gradiente Coniugato");

    System.out.print("Scelta: ");
    int choice = sc.nextInt();

    // **************************
    // CARICAMENTO MATRICE E VETTORE 
    // **************************
    DMatrixSparseCSC A = MatrixMarketReader.readMatrixMarketFile(path);
    int n = A.numRows;

    DMatrixRMaj xTrue = new DMatrixRMaj(n, 1);
    for (int i = 0; i < n; i++) {
        xTrue.set(i, 0, 1.0);
    }

    DMatrixRMaj b = new DMatrixRMaj(n, 1);
    org.ejml.sparse.csc.CommonOps_DSCC.mult(A, xTrue, b);

    // **************************
    // SELEZIONE DEL SOLVER
    // **************************
    IterativeSolver solver = switch (choice) {
        case 1 -> new JacobiSolver();
        case 2 -> new GaussSeidel();
        case 3 -> new Gradient();
        case 4 -> new ConjugateGradient();
        default -> throw new IllegalArgumentException("Metodo non valido");
    };

    // **************************
    // ESECUZIONE
    // **************************
    SolverResult result = solver.solve(A, b, tol, maxIter);

    double error =
            UtilsOperations.relativeError(xTrue, result.getSolution());

    // **************************
    // OUTPUT
    // **************************
    System.out.println("\n========== RISULTATO ==========");
    System.out.println("Matrice: " + matrices[m - 1]);
    System.out.println("Tolleranza: " + tol);
    System.out.println("Metodo: " + choice);
    System.out.println("--------------------------------");
    System.out.println("Iterazioni: " + result.getIterations());
    System.out.println("Residuo relativo: " + result.getRelativeResidual());
    System.out.println("Errore relativo: " + error);
    System.out.println("Tempo: " + result.getExecutionTime() + " s");
    System.out.println("Convergenza: " + result.isConverged());
    System.out.println("================================");
    }
}