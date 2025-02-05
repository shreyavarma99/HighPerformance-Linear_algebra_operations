#include "assignment1.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

int test_ger(int nrepeats, int first, int last, int inc)
{
    int size, irep;
    int m, n;
    int csA, rsA;

    double *A_my, *A_ref;
    double *x, *y;

    double t_ref = DBL_MAX;
    double t = DBL_MAX;
    double t_start;

    double alpha = 1.0;
    double maxdiff = 0.0;

    printf("%% --------- DGER --------- \n");

    for (size = last; size >= first; size -= inc){
        m = size;
        n = size / 2;
        csA = n;
        rsA = 1;

        int incx = size;
        int incy = size;

        A_my = (double *)malloc(m * n * sizeof(double));
        A_ref = (double *)malloc(m * n * sizeof(double));
        x = (double *)malloc(m * incx * sizeof(double));
        y = (double *)malloc(n * incy * sizeof(double));

        bli_drandm(0, BLIS_DENSE, m, n, A_my, rsA, csA);
        memcpy(A_ref, A_my, m * n * sizeof(double));
        bli_drandv(m, x, incx);
        bli_drandv(n, y, incy);

        for (irep = 0; irep < nrepeats; irep++){
            t_start = bli_clock();
            bli_dger(
                BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE,
                m, n,
                &alpha,
                x, incx,
                y, incy,
                A_ref, rsA, csA);
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }

        double gflops_ref = (2.0 * m * n) / (t_ref * 1.0e9);

        for (irep = 0; irep < nrepeats; irep++){
            t_start = bli_clock();
            shpc_dger(m, n, x, incx, y, incy, A_my, rsA, csA);
            t = bli_clock_min_diff(t, t_start);
        }

        double gflops = (2.0 * m * n) / (t * 1.0e9);

        double diff = shpc_maxabsdiff(m, n, A_my, rsA, csA, A_ref, rsA, csA);
        maxdiff = max(diff, maxdiff);

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
