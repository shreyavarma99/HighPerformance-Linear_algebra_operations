#include "assignment1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "blis.h"


int test_axpy(int nrepeats, int first, int last, int inc) {
    int size, irep;
    double *x, *y_my, *y_ref;
    double alpha = 2.0;
    double t_ref = DBL_MAX, t = DBL_MAX, t_start;
    double diff, maxdiff = 0.0;
    const double tolerance = 1e-12;
    double gflops_ref, gflops;


    printf("%% --------- AXPY --------- \n");
    printf("data_axpy");
    printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
           (unsigned long)(last - first) / inc + 1,
           (unsigned long)0, 0.0, 0.0, 0.0);


    for (size = last; size >= first; size -= inc) {
        x = (double *)malloc(size * sizeof(double));
        y_my = (double *)malloc(size * sizeof(double));
        y_ref = (double *)malloc(size * sizeof(double));


        for (int i = 0; i < size; i++) {
            x[i] = (double)(rand() % 100) / 10.0;
            y_my[i] = y_ref[i] = (double)(rand() % 100) / 10.0;
        }


        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            bli_daxpyv(BLIS_NO_CONJUGATE, size, &alpha, x, 1, y_ref, 1);
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }

        gflops_ref = (2.0 * size) / (t_ref * 1.0e9);


        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            shpc_daxpy(size, &alpha, x, 1, y_my, 1);
            t = bli_clock_min_diff(t, t_start);
        }

        gflops = (2.0 * (double) size) / (t * 1.0e9);

        maxdiff = 0.0;
        for (int i = 0; i < size; i++) {
            diff = fabs(y_my[i] - y_ref[i]);
            if (diff > maxdiff) {
                maxdiff = diff;
            }
        }


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
