#include "assignment1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "blis.h"


// Function to test the AXPY implementation
int test_axpy(int nrepeats, int first, int last, int inc) {
    int size, irep;
    double *x, *y_my, *y_ref;
    double alpha = 2.0;
    double t_ref = DBL_MAX, t = DBL_MAX, t_start;
    double diff, maxdiff = 0.0;
    const double tolerance = 1e-12;


    printf("%% --------- AXPY --------- \n");
    printf("data_axpy");
    printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
           (unsigned long)(last - first) / inc + 1,
           (unsigned long)0, 0.0, 0.0, 0.0);


    // Loop through various sizes of vectors
    for (size = last; size >= first; size -= inc) {
        // Allocate vectors
        x = (double *)malloc(size * sizeof(double));
        y_my = (double *)malloc(size * sizeof(double));
        y_ref = (double *)malloc(size * sizeof(double));


        // Initialize vectors with random values
        for (int i = 0; i < size; i++) {
            x[i] = (double)(rand() % 100) / 10.0;
            y_my[i] = y_ref[i] = (double)(rand() % 100) / 10.0;
        }


        // double sum_check = 0.0;
        // for (int i = 0; i < size; i++) {
        //     printf("%8.4f + ", y_ref[i]);
        //     sum_check += y_ref[i];
        // }


        // printf("  >> Sum of y before: %f\n", sum_check);


        //  sum_check = 0.0;
        // for (int i = 0; i < size; i++) {
        //     printf("%8.4f + ", x[i]);
        //     sum_check += x[i];
        // }


        // printf("  >> Sum of x: %f\n", sum_check);


        // Run BLIS reference implementation
        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            bli_daxpyv(BLIS_NO_CONJUGATE, size, &alpha, x, 1, y_ref, 1);
            t_ref = bli_clock_min_diff(t_ref, t_start);
        }


        //  sum_check = 0.0;
        // for (int i = 0; i < size; i++) {
        //     sum_check += y_ref[i];
        // }
        // printf("Sum of y my: %f\n", sum_check);


        // Run my implementation
        for (irep = 0; irep < nrepeats; irep++) {
            t_start = bli_clock();
            shpc_daxpy(size, &alpha, x, 1, y_my, 1);
            t = bli_clock_min_diff(t, t_start);
        }
        // sum_check = 0.0;
        // for (int i = 0; i < size; i++) {
        //     sum_check += y_my[i];
        // }
        // printf("Sum of y ref: %f\n", sum_check);


        // Compare the results
        maxdiff = 0.0;
        for (int i = 0; i < size; i++) {
            diff = fabs(y_my[i] - y_ref[i]);
            if (diff > maxdiff) {
                maxdiff = diff;
            }
        }


        // Print the results
        printf("data_axpy");
        printf("( %4lu, 1:5 ) = [ %5lu %8.2f %8.2f %15.4e ];\n",
               (unsigned long)(size - first) / inc + 1,
               (unsigned long)size, t_ref, t, maxdiff);




        // Free allocated memory
        free(x);
        free(y_my);
        free(y_ref);
    }


    return 0;
}
