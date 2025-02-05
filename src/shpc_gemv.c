void shpc_dgemv(int m, int n, double *A, int rsA, int csA, double *x, int incx, double *y, int incy) {
    for (int i = 0; i < m; i++) {
        int y_index = i * incy;
        double temp = 0.0;

        // Accumulate the row-wise dot product for y[i]
        for (int j = 0; j < n; j++) {
            int A_index = i * rsA + j * csA;
            int x_index = j * incx;
            temp += A[A_index] * x[x_index];
        }

        // Update y[i]
        y[y_index] += temp;
    }
}
