#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>

int main() {
    int N = 10;
    double a = 3.0;
    double x_value = 0.1;
    double y_value = 7.1;

    // allocate vectors and initialize
    gsl_vector *x = gsl_vector_alloc(N);
    gsl_vector *y = gsl_vector_alloc(N);
    for (int i = 0; i < N; i++) {
        gsl_vector_set(x, i, x_value);
        gsl_vector_set(y, i, y_value);
    }

    gsl_vector_axpby(a, x, 1.0, y);
    printf("[info] look at first result: d[0] = %f\n", gsl_vector_get(y, 0));

    // dump to file ... (skipping)

    // free memory
    gsl_vector_free(x);
    gsl_vector_free(y);

    return EXIT_SUCCESS;
}
