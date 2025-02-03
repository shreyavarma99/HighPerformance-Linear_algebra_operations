#include "assignment1.h"

void shpc_dgemm( int m, int n, int k,                                            
                    double *A, int rsA, int csA,                                
                    double *B, int rsB, int csB,                                
                    double *C, int rsC, int csC )
{

    // Your code goes here. 
    for(int j = 0; j < n; j++) { //columns
        for(int i = 0; i < m; i++) { //rows
            double sum = 0.0;
            for(int l = 0; l < k; l++){
               
                int A_index = i * rsA + l * csA;
                int B_index = l * rsB + j * csB;
                sum += A[A_index] * B[B_index];
            }
            int C_index = i * rsC + j * csC;
            C[C_index] += sum;
        }
    }

}   
