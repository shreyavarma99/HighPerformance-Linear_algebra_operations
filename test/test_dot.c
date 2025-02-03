#include "assignment1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "blis.h"


int test_dot(int nrepeats, int first, int last, int inc) {
    int size, irep;
    double *x, *y;
    double rho_ref = 0.0, rho = 0.0;
    double t_ref = DBL_MAX, t = DBL_MAX, t_start;
    double diff, maxdiff = 0.0;
    const double tolerance = 1e-12;


    printf("%% --------- DDOT --------- \n");
    printf("data_ddot");
    printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
           (unsigned long)(last - first) / inc + 1,
           (unsigned long)0, 0.0, 0.0, 0.0);


    for (size = last; size >= first; size -= inc) {
        // Allocate vectors
        x = (double *)malloc(size * sizeof(double));
        y = (double *)malloc(size * sizeof(double));


        for (int i = 0; i < size; i++) {
            x[i] = (double)(rand() % 100) / 10.0;
            y[i] = (double)(rand() % 100) / 10.0;
        }


        for (irep = 0; irep < nrepeats; irep++) {
            rho_ref = 0.0;
            t_start = bli_clock();
            bli_ddotv(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, size, x, 1, y, 1, &rho_ref);
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }


        for (irep = 0; irep < nrepeats; irep++) {
            rho = 0.0;
            t_start = bli_clock();
            shpc_ddot(size, x, 1, y, 1, &rho);
            t = bli_clock_min_diff(t, t_start);
        }


        diff = fabs(rho - rho_ref);
        maxdiff = fmax(diff, maxdiff);


        printf("data_ddot");
        printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
               (unsigned long)(size - first) / inc + 1,
               (unsigned long)size, t_ref, t, diff);

        free(x);
        free(y);
    }


    return 0;
}
