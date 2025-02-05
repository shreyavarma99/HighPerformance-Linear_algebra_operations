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

        bli_drandm(0, BLIS_DENSE, m, n, A_my, rsA, csA);
        memcpy(A_ref, A_my, m * n * sizeof(double));
        bli_drandv(m, x, 1);
        bli_drandv(n, y, 1);

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

        double gflops_ref = (2.0 * m * n) / (t_ref * 1.0e9);


        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            shpc_dger(m, n, x, 1, y, 1, A_my, rsA, csA);
            t = bli_clock_min_diff(t, t_start);
        }

        double gflops = (2.0 * m * n) / (t * 1.0e9);


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
               (unsigned long)m, (unsigned long)n, (unsigned long)0, gflops_ref, gflops, maxdiff);

        free(A_my);
        free(A_ref);
        free(x);
        free(y);
    }

    return 0;
}
