#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// example constants
#define EXAMPLE_ARRAY_SIZE 100000000 // also test 10, 1000000 
#define EXAMPLE_EXPECTED_RESULT 7.4
#define DOUBLE_EPSILON 1e-14

void daxpy_not_inplace(int n, double a, double *px, double *py, double *pd) {
    // computes pd = a*px + py
    // where int n is the size of arrays px, py, pd
    int i, m;

    // return on trivial cases
    if (n<=0) return;
    if (a==0.0) {
        for (i=0; i<n; i++) pd[i] = 0.0;
    }

    // clean-up loop
    m = n%4;
    if (m!=0) {
        for (i=0; i<m; i++) {
            pd[i] = py[i] + a * px[i];
        }
        if (n<4) return;
    }
    // loop unrolling in units of 4
    // (assuming 256-bit registers)
    for (i=m; i<n; i+=4) {
        pd[i]   = py[i] + a*px[i];
        pd[i+1] = py[i+1] + a*px[i+1];
        pd[i+2] = py[i+2] + a*px[i+2];
        pd[i+3] = py[i+3] + a*px[i+3];
    }

}


int main() {
    double a;
    double *px, *py, *pd;
    int i = 0;
    bool test = true;

    // allocating arrays
    printf("allocating...\n");
    px = (double*)malloc(EXAMPLE_ARRAY_SIZE * sizeof(double));
    py = (double*)malloc(EXAMPLE_ARRAY_SIZE * sizeof(double));
    pd = (double*)malloc(EXAMPLE_ARRAY_SIZE * sizeof(double));

    if (px == NULL || py == NULL || pd == NULL) {
        printf("[error] memory allocation failed\n");
        return 1;
    }

    // initialize arrays & constants for this test
    for (i = 0; i < EXAMPLE_ARRAY_SIZE; i++) {
        px[i] = 0.1;
        py[i] = 7.1;
    }
    a = 3;

    // daxpy (note: output on pd)
    printf("computing...\n");
    daxpy_not_inplace(EXAMPLE_ARRAY_SIZE, a, px, py, pd);

    // testing the result
    printf("checking result...");
    for (i = 0; i < EXAMPLE_ARRAY_SIZE; i++)
        if (fabs(pd[i] - EXAMPLE_EXPECTED_RESULT) >= DOUBLE_EPSILON) break;

    // printing the result
    if (i < EXAMPLE_ARRAY_SIZE) {
        printf(" test failed!\n");
        printf("  expected result is %f, but for instance d[%d]=%f\n", EXAMPLE_EXPECTED_RESULT, i, pd[i]);
    }
    else printf(" all elements are correct\n");

    // free allocated memory
    free(px);
    free(py);
    free(pd);

    return 0;
}
