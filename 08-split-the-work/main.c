#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

double* daxpy_explicit_chunked(double *d, double a, double *x, double *y, int n, int chunk_size) {
    // trivial case, a=0 so inplace update is a no-op
    if (n <= 0 || a == 0.0) return NULL;

    int current_start, current_end, chunk_index, i;

    int num_chunks = (int)ceil( (double)n / (double)chunk_size );
    double* partial_chunk_sum = malloc( num_chunks * sizeof(double) );

    for (chunk_index = 0; chunk_index < num_chunks; chunk_index++) {
        current_start = chunk_index * chunk_size;
        current_end = current_start + chunk_size;
        if (current_end > n) current_end = n;

        for (i = current_start; i < current_end; i++) {
            d[i] += a * x[i] + y[i];
            partial_chunk_sum[chunk_index] += d[i];
        }
    }
    return partial_chunk_sum;
}



int main() {
    int n = 1000;
    int chunk_size = 8;
    double a = 2.0;

    double x[n], y[n], d[n];
    double *partial_chunk_sum, *d_chunk;

    srand(12345); // fix the seed for reproducibility

    // fill x and y with random numbers in [-1,1]
    for (int i = 0; i < n; i++) {
        x[i] = ((double) rand() / RAND_MAX) * 2.0 - 1.0;
        y[i] = ((double) rand() / RAND_MAX) * 2.0 - 1.0;
    }

    // DAXPY, but not inplace to preserve y
    double sum_d = 0.0;
    for (int i = 0; i < n; i++) {
        d[i] = a * x[i] + y[i];
        sum_d += d[i];
    }

    // chunked DAXPY (allocates partial_chunk_sum internally)
    d_chunk = malloc(n * sizeof(double));
    partial_chunk_sum = daxpy_explicit_chunked(d_chunk, a, x, y, n, chunk_size);

    // 1) check that d is the same as the original code
    for (int i = 0; i < n; i++)
        assert(fabs(d[i] - d_chunk[i]) < 1e-12);
    printf(" [test] daxpy chunked vs original passed\n");

    // 2) compute total sum from partial sums
    double sum_dchunked = 0.0;
    for (int i = 0; i < (int)ceil((double)n / chunk_size); i++)
        sum_dchunked += partial_chunk_sum[i];
    assert(fabs(sum_d - sum_dchunked) < 1e-10);
    printf(" [test] sum vs sum chunked passed: %f vs %f\n", sum_d, sum_dchunked);

    // cleanup
    free(partial_chunk_sum);
    free(d_chunk);
    return 0;
}

