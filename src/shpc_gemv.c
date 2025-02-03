void shpc_dgemv(int m, int n, double *A, int rsA, int csA, double *x, int incx, double *y, int incy) {
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            int A_index = i * rsA + j * csA;
            int x_index = j * incx;
            sum += A[A_index] * x[x_index];
        }
        int y_index = i * incy;
        y[y_index] += sum;
    }
}
