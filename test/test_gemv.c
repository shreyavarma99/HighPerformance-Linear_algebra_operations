#include "assignment1.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

int test_gemv(int nrepeats, int first, int last, int inc) {
    int size, irep;
    int m, n;
    int csA, rsA;

    double *A, *x, *y_my, *y_ref;

    double t_ref = DBL_MAX;
    double t = DBL_MAX;
    double t_start;

    double gflops_ref, gflops;
    double alpha = 1.0, beta = 1.0;
    int incx, incy;

    printf("%% --------- DGEMV --------- \n");

    for (size = last; size >= first; size -= inc) {
        m = size;
        n = size / 2;
        csA = n;
        rsA = 1;

        //smaller stride
        incx = size;
        incy = size;

        A = (double *)malloc(m * n * sizeof(double));
        x = (double *)malloc(n * incx * sizeof(double));
        y_my = (double *)malloc(m * incy * sizeof(double));
        y_ref = (double *)malloc(m * incy * sizeof(double));

        bli_drandm(0, BLIS_DENSE, m, n, A, rsA, csA);
        bli_drandv(n, x, incx);
        bli_drandv(m, y_my, incy);
        memcpy(y_ref, y_my, m * incy * sizeof(double));

        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            bli_dgemv(
                BLIS_NO_TRANSPOSE,
                BLIS_NO_CONJUGATE,
                m, n,
                &alpha,
                A, rsA, csA,
                x, incx,
                &beta,
                y_ref, incy
            );
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }

        gflops_ref = (2.0 * m * n) / (t_ref * 1.0e9);

        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            shpc_dgemv(m, n, A, rsA, csA, x, incx, y_my, incy);
            t = bli_clock_min_diff(t, t_start);
        }

        gflops = (2.0 * m * n) / (t * 1.0e9);

        double maxdiff = shpc_maxabsdiff(m, 1, y_my, incy, 1, y_ref, incy, 1);

        printf("data_gemv");
        printf("( %4lu, 1:6 ) = [ %5lu %5lu %5lu %8.2f %8.2f %15.4e ];\n",
               (unsigned long)(size - first) / inc + 1,
               (unsigned long)m, (unsigned long)n, (unsigned long)0, gflops_ref, gflops, maxdiff);

        free(A);
        free(x);
        free(y_my);
        free(y_ref);
    }

    return 0;
}
