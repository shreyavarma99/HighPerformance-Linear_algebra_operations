#include "assignment1.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

int test_ger(int nrepeats, int first, int last, int inc) {
    int size, irep;
    int m, n;
    int csA, rsA;

    double *A_my, *A_ref;
    double *x, *y;

    double t_ref = DBL_MAX;
    double t = DBL_MAX;
    double t_start;

    double diff, maxdiff = 0.0;
    const double tolerance = 1e-12;

    double alpha = 1.0;

    printf("%% --------- DGER --------- \n");
    printf("data_ger");
    printf("( %4lu, 1:6 ) = [ %5lu %5lu %5lu %8.2f %8.2f %15.4e ];\n",
           (unsigned long)(last - first) / inc + 1,
           (unsigned long)0, (unsigned long)0, (unsigned long)0, 0.0, 0.0, 0.0);

    for (size = last; size >= first; size -= inc) {
        m = size;
        n = size / 2;
        csA = n;
        rsA = 1;

        A_my = (double *)malloc(m * n * sizeof(double));
        A_ref = (double *)malloc(m * n * sizeof(double));
        x = (double *)malloc(m * sizeof(double));
        y = (double *)malloc(n * sizeof(double));

        for (int i = 0; i < m * n; i++) {
            A_my[i] = A_ref[i] = (double)(rand() % 100) / 10.0;
        }

        for (int i = 0; i < m; i++) {
            x[i] = (double)(rand() % 100) / 10.0;
        }

        for (int i = 0; i < n; i++) {
            y[i] = (double)(rand() % 100) / 10.0;
        }

        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            bli_dger(
                BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE,
                m, n,
                &alpha,
                x, 1, 
                y, 1,
                A_ref, rsA, csA
            );
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }

        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            shpc_dger(m, n, x, 1, y, 1, A_my, rsA, csA);
            t = bli_clock_min_diff(t, t_start);
        }

        maxdiff = 0.0;
        for (int i = 0; i < m * n; i++) {
            diff = fabs(A_my[i] - A_ref[i]);
            if (diff > maxdiff) {
                maxdiff = diff;
            }
        }

        printf("data_ger");
        printf("( %4lu, 1:6 ) = [ %5lu %5lu %5lu %8.2f %8.2f %15.4e ];\n",
               (unsigned long)(size - first) / inc + 1,
               (unsigned long)m, (unsigned long)n, (unsigned long)0, t_ref, t, maxdiff);

        free(A_my);
        free(A_ref);
        free(x);
        free(y);
    }

    return 0;
}
