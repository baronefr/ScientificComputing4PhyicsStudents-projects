#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

// target routine to test: DAXPY (Double-precision AX Plus Y)
void daxpy(double a, double *x, double *y, int n) {
    // NOTE: result is inplace in y
    for (int i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}

// Gaussian random number generator
double gaussian_random() {
    double u1 = ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
    double u2 = ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}



// ---------- Tests ----------

void test_daxpy_small() {
    double x[3] = {1.0, 2.0, 3.0};
    double y[3] = {4.0, 5.0, 6.0};
    const double tolerance = 1e-15;

    daxpy(2.0, x, y, 3); // expected y = {6.0, 9.0, 12.0}
    // note: inplace operation on y

    assert(fabs(y[0] - 6.0)  < tolerance);
    assert(fabs(y[1] - 9.0)  < tolerance);
    assert(fabs(y[2] - 12.0) < tolerance);
    printf("[test] daxpy_small passed\n");
}


void test_daxpy_zero() {
    // test with a = 0, should return y
    double x[3] = {1.0, 2.0, 3.0};
    double y[3] = {4.0, 5.0, 6.0};
    const double tolerance = 1e-15;

    daxpy(0.0, x, y, 3); // expected y = {6.0, 9.0, 12.0}
    // note: inplace operation on y

    assert(fabs(y[0] - 4.0)  < tolerance);
    assert(fabs(y[1] - 5.0)  < tolerance);
    assert(fabs(y[2] - 6.0)  < tolerance);
    printf("[test] daxpy_zero passed\n");
}


void test_daxpy_rand() {
    // test with random data
    int test_size = (int)1e5;
    double *x = malloc(test_size * sizeof(double));
    double *y = malloc(test_size * sizeof(double));

    // fill with random data
    for (int i = 0; i < test_size; i++) {
        x[i] = gaussian_random();
        y[i] = gaussian_random();
    }
    double a = gaussian_random();

    // compute expected result
    double *y_expected = malloc(test_size * sizeof(double));
    for (int i = 0; i < test_size; i++) {
        y_expected[i] = y[i] + a * x[i];
    }
    // call daxpy
    daxpy(a, x, y, test_size);

    // compare
    const double tolerance = 1e-15;
    for (int i = 0; i < test_size; i++) {
        assert(fabs(y[i] - y_expected[i]) < tolerance);
    }
    printf("[test] daxpy_rand passed\n");

    free(x);
    free(y);
    free(y_expected);
}


int main() {
    srand(time(NULL));
    test_daxpy_small();
    test_daxpy_zero();
    test_daxpy_rand();
    printf(" all tests passed!\n");
    return 0;
}
