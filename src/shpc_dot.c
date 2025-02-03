#include <stdio.h>
void shpc_ddot (
int n ,
double *x , int incx ,
double *y , int incy ,
double * rho
){
    double sum = 0;
    for(int i = 0; i < n; i++){
        sum += x[i * incx] * y[i * incy];
    }
    *rho = sum;
}
