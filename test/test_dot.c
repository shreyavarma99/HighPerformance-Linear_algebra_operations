#include "assignment1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "blis.h"

int test_dot(int nrepeats, int first, int last, int inc)
{
    int size, irep;
    double *x, *y;
    int incX, incY;

    double alpha = 1.0, beta = 1.0;

    double done = 1.0;
    double t_ref = DBL_MAX;
    double t = DBL_MAX;
    double t_start;

    double gflops_ref, gflops;

    double diff, maxdiff = 0.0;

    printf("%% --------- DDOT --------- \n");
    printf("data_ddot");
    printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
           (unsigned long)(last - first) / inc + 1,
           (unsigned long)0, 0.0, 0.0, 0.0);

    for (size = last; size >= first; size -= inc)
    {
        incX = size;
        incY = size;

        double rho = (double)(rand() % 100) / 10.0;
        double rho_my = rho;
        double rho_ref = rho;

        // Allocate vectors
        x = (double *)malloc(size * incX * sizeof(double));
        y = (double *)malloc(size * incY * sizeof(double));

        bli_drandv(size, x, incX);
        bli_drandv(size, y, incY);

        for (irep = 0; irep < nrepeats; irep++)
        {
            t_start = bli_clock();
            bli_ddotxv(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, size, &alpha, x, incX, y, incY, &beta, &rho_ref);
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }

        gflops_ref = (2.0 * size) / (t_ref * 1.0e9);

        for (irep = 0; irep < nrepeats; irep++)
        {
            t_start = bli_clock();
            shpc_ddot(size, x, incX, y, incY, &rho_my);
            t = bli_clock_min_diff(t, t_start);
        }

        gflops = (2.0 * (double)size) / (t * 1.0e9);
        diff = dabs(rho_my - rho_ref);
        maxdiff = max(diff, maxdiff);

        printf("data_ddot");
        printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
               (unsigned long)(size - first) / inc + 1,
               (unsigned long)size, gflops_ref, gflops, maxdiff);

        free(x);
        free(y);
    }

    return 0;
}
