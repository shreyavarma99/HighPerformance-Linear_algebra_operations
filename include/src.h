
/* Computes C = AB + C */
void shpc_dgemm( int m, int n, int k, 
                    double *A, int rsA, int csA, 
                    double *B, int rsB, int csB, 
                    double *C, int rsC, int csC );
void shpc_daxpy (int n , double * alpha , double *x , int incx , double *y , int incy);
void shpc_ddot (int n , double *x , int incx , double *y , int incy , double * rho);
void shpc_dgemv(int m, int n, double *A, int rsA, int csA, double *x, int incx, double *y, int incy);
void shpc_dger (int m , int n , double *x , int incx , double *y , int incy , double *A , int rsA , int csA);


                
