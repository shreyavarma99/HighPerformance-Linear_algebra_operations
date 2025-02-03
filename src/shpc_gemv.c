//i on the inside made error go lower

void shpc_dgemv(int m, int n, double *A, int rsA, int csA, double *x, int incx, double *y, int incy) {
    for (int j = 0; j < n; j++) {  // Loop over columns first
        int x_index = j * incx;
        for (int i = 0; i < m; i++) {  // Loop over rows inside
            int A_index = i * rsA + j * csA;
            int y_index = i * incy;
            y[y_index] += A[A_index] * x[x_index];
        }
    }
}