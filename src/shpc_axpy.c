void shpc_daxpy (int n , double * alpha , double *x , int incx , double *y , int incy){
    for(int i = 0; i < n; i++){
        int x_index = i * incx;
        int y_index = i * incy;


        y[y_index] += (x[x_index] * (*alpha));
    }
}
