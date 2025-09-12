
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// example constants
#define EXAMPLE_ARRAY_SIZE 1000000
#define EXAMPLE_A 3.14


void daxpy(int n, double a, double *px, double *py) {
    // computes a*px + py, in-place on py
    // where int n is the size of arrays px, py
    int i, m;

    // return on trivial cases
    if (n<=0) return;
    if (a==0.0) return;

    // clean-up loop
    m = n%4;
    if (m!=0) {
        for (i=0; i<m; i++) {
            py[i] += a * px[i];
        }
        if (n<4) return;
    }
    // loop unrolling in units of 4
    // (assuming 256-bit registers)
    for (i=m; i<n; i+=4) {
        py[i]   += a*px[i];
        py[i+1] += a*px[i+1];
        py[i+2] += a*px[i+2];
        py[i+3] += a*px[i+3];
    }

}


int main() {
    double a = EXAMPLE_A;
    double *px, *py;

    // rng seed
    srand(time(NULL));

    // allocating arrays
    px = (double*)malloc(EXAMPLE_ARRAY_SIZE * sizeof(double));
    py = (double*)malloc(EXAMPLE_ARRAY_SIZE * sizeof(double));

    if (px == NULL || py == NULL) {
        printf("[error] memory allocation failed\n");
        return 1;
    }

    // initialize arrays with random doubles in [0, 1)
    for (int i = 0; i < EXAMPLE_ARRAY_SIZE; i++) {
        px[i] = (double)rand() / RAND_MAX;
        py[i] = (double)rand() / RAND_MAX;
    }

    // daxpy (note: in-place on py)
    daxpy(EXAMPLE_ARRAY_SIZE, a, px, py);

    // print something
    printf("Printing first 10 elements:\n");
    for (int i = 0; i < 10; i++) {
        printf(" y[%d] = %f\n", i, py[i]);
    }

    // free allocated memory
    free(px);
    free(py);

    return 0;
}
