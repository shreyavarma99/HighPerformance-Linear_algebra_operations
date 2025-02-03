void shpc_dger (int m , int n , double *x , int incx , double *y , int incy , double *A , int rsA , int csA){
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            int x_index = incx * i;
            int y_index = incy * j;
            A[rsA * i + j * csA] += x[x_index] * y[y_index];
        }
    }
}