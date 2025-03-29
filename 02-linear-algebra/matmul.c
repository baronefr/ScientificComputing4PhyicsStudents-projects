#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void matmul(double *A, double *B, double *C, int n) {
    // naive way to do matmul ...
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            C[i*n + j] = 0.0;
            for (int k = 0; k < n; k++)
                C[i*n + j] += A[i*n + k] * B[k*n + j];
        }
}

bool allclose(double *C, double expected, int n, double epsilon) {
    for (int i = 0; i < n*n; i++)
        if (fabs(C[i] - expected) >= epsilon) return false;
    return true;
}

int main() {
    int N[] = {10, 100, 10000};
    int num_tests = sizeof(N) / sizeof(N[0]);
    double *A, *B, *C;
    int nn;
    bool test;

    for (int t = 0; t < num_tests; t++) {
        nn = N[t];
        printf("testing for N=%d ...\n", nn);

        A = (double*)malloc(nn * nn * sizeof(double));
        B = (double*)malloc(nn * nn * sizeof(double));
        C = (double*)malloc(nn * nn * sizeof(double));

        if (A == NULL || B == NULL || C == NULL) {
            printf("[error] memory allocation failed\n");
            return 1;
        }

        for (int i = 0; i < nn * nn; i++) {
            A[i] = 0.3;
            B[i] = 7.1;
        }

        matmul(A, B, C, nn);

        printf(" checking the result... ");
        test = allclose(C, 0.3 * 7.1 * nn, nn, 1e-8);
        if (test) printf(" passed\n");
        else {
            printf(" failed\n");
            printf("  C[0] = %.15f but expected %.15f\n", C[0], 0.3 * 7.1 * nn);
            return 1;
        }

        free(A);
        free(B);
        free(C);
    }

    return 0;
}
