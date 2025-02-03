#include "assignment1.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

// Function to test GEMV
int test_gemv(int nrepeats, int first, int last, int inc) {
    int size, irep;
    int m, n;
    int csA, rsA;

    double *A, *x, *y_my, *y_ref;

    double t_ref = DBL_MAX;
    double t = DBL_MAX;
    double t_start;

    double diff, maxdiff = 0.0;
    const double tolerance = 1e-12;

    double alpha = 1.0, beta = 1.0;

    printf("%% --------- DGEMV --------- \n");
    printf("data_gemv");
    printf("( %4lu, 1:6 ) = [ %5lu %5lu %5lu %8.2f %8.2f %15.4e ];\n",
           (unsigned long)(last - first) / inc + 1,
           (unsigned long)0, (unsigned long)0, (unsigned long)0, 0.0, 0.0, 0.0);

    // Loop through different sizes
    for (size = last; size >= first; size -= inc) {
        m = size;            // Rows of A
        n = size / 2;        // Columns of A (or adjust as needed)
        csA = n;             // Column stride (elements per row)
        rsA = 1;             // Row stride

        // Allocate matrix A and vectors x, y_my, y_ref
        A = (double *)malloc(m * n * sizeof(double));
        x = (double *)malloc(n * sizeof(double));
        y_my = (double *)malloc(m * sizeof(double));
        y_ref = (double *)malloc(m * sizeof(double));

        // Initialize matrix and vectors with random values
        for (int i = 0; i < m * n; i++) {
            A[i] = (double)(rand() % 100) / 10.0;
        }

        for (int i = 0; i < n; i++) {
            x[i] = (double)(rand() % 100) / 10.0;
        }

        for (int i = 0; i < m; i++) {
            y_my[i] = y_ref[i] = 0.0;
        }

        // Run BLIS reference implementation
        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            bli_dgemv(
                BLIS_NO_TRANSPOSE,               // Matrix A is not transposed
                BLIS_NO_CONJUGATE,               // No conjugation for x
                m, n,                            // Matrix dimensions
                &alpha,                          // Alpha scalar
                A, rsA, csA,                     // Matrix A with row/column strides
                x, 1,                            // Vector x with stride 1
                &beta,                           // Beta scalar
                y_ref, 1                         // Output vector y_ref with stride 1
            );
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }

        // Run your implementation
        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            shpc_dgemv(m, n, A, rsA, csA, x, 1, y_my, 1);
            t = bli_clock_min_diff(t, t_start);
        }

        // Compare results
        maxdiff = 0.0;
        for (int i = 0; i < m; i++) {
            diff = fabs(y_my[i] - y_ref[i]);
            if (diff > maxdiff) {
                maxdiff = diff;
            }
        }

        // Print results
        printf("data_gemv");
        printf("( %4lu, 1:6 ) = [ %5lu %5lu %5lu %8.2f %8.2f %15.4e ];\n",
               (unsigned long)(size - first) / inc + 1,
               (unsigned long)m, (unsigned long)n, (unsigned long)0, t_ref, t, maxdiff);

        // Free allocated memory
        free(A);
        free(x);
        free(y_my);
        free(y_ref);
    }

    return 0;
}
