#include "assignment1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "blis.h"

int test_axpy(int nrepeats, int first, int last, int inc) {
    int size, irep;
    double *x, *y_my, *y_ref;
    double t_ref = DBL_MAX, t = DBL_MAX, t_start;
    double maxdiff = 0.0;
    double gflops_ref, gflops;
    int incx, incy;

    printf("%% --------- AXPY --------- \n");
    printf("data_axpy");
    printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
           (unsigned long)(last - first) / inc + 1,
           (unsigned long)0, 0.0, 0.0, 0.0);

    for (size = last; size >= first; size -= inc) {
        double alpha = (double)(rand() % 100) / 10;

        incx = size;
        incy = size;

        x = (double *)malloc(incx * size * sizeof(double));
        y_my = (double *)malloc(incy * size * sizeof(double));
        y_ref = (double *)malloc(incy * size * sizeof(double));

        bli_drandv(size, x, incx);
        bli_drandv(size, y_my, incy);
        memcpy(y_ref, y_my, size * incy * sizeof(double));

        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            bli_daxpyv(BLIS_NO_CONJUGATE, size, &alpha, x, incx, y_ref, incy);
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }

        gflops_ref = (2.0 * size) / (t_ref * 1.0e9);

        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            shpc_daxpy(size, &alpha, x, incx, y_my, incy);
            t = bli_clock_min_diff(t, t_start);
        }

        gflops = (2.0 * size) / (t * 1.0e9);

        double diff = shpc_maxabsdiff(size, 1, y_my, incy, 1, y_ref, incy, 1);
        maxdiff = max(diff, maxdiff);

        printf("data_axpy");
        printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
               (unsigned long)(size - first) / inc + 1,
               (unsigned long)size, gflops_ref, gflops, maxdiff);

        free(x);
        free(y_my);
        free(y_ref);
    }

    return 0;
}
