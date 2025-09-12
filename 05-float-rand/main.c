#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <assert.h>

double forloop_sum(double *vec, int n) {
    // trivially sum the elements in a loop
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += vec[i];
    }
    return sum;
}

double kahan_sum(double *vec, int n) {
    // Kahan summation algorithm
    double sum = 0.0;
    double c = 0.0; // c for compensation
    double y = 0.0;
    double t = 0.0;

    for (int i = 0; i < n; i++) {
        y = vec[i] - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
        //printf("i=%d, y=%.17g, t=%.17g, c=%.17g, sum=%.17g\n", i, y, t, c, sum);
    }
    return sum;
}

// Kahan–Babuška–Neumaier summation
double kbn_sum(const double *a, int n) {
    double sum = a[0];
    double c = 0.0; // correction

    for (int i = 1; i < n; i++) {
        double t = sum + a[i];
        if (fabs(sum) >= fabs(a[i])) {
            c += (sum - t) + a[i];
        } else {
            c += (a[i] - t) + sum;
        }
        sum = t;
    }
    return sum + c;
}

double gsl_wrapped_sum(double *vec, int n) {
    // first I create a gsl_vector_view from the array
    gsl_vector_view v = gsl_vector_view_array(vec, n);
    // then I call gsl_vector_sum
    return gsl_vector_sum(&v.vector);
}


// ---------- Part (b) ----------
// DAXPY: y = a*x + y
void daxpy(double a, double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}

// Gaussian RNG (Box-Muller)
double gaussian_random() {
    double u1 = ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
    double u2 = ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// create a vector of n Gaussian random numbers
void fill_gaussian(double *vec, int n) {
    for (int i = 0; i < n; i++) {
        vec[i] = gaussian_random();
    }
}

// Test: check sum of (x+y) by direct calculation
double test_sum(double *x, double *y, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        s += x[i] + y[i];
    }
    return s;
}

// ---------- Statistics ----------
double compute_mean(const double *vec, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += vec[i];
    return sum / n;
}

double compute_std(const double *vec, int n, double mean) {
    double sumsq = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = vec[i] - mean;
        sumsq += diff * diff;
    }
    return sqrt(sumsq / n);
}

int main() {
    // subtask a -----------------------
    double vec[4] = {1.0, 1.0e16, -1.0e16, -0.5};

    double sum_for = forloop_sum(vec, 4);
    double sum_gsl = gsl_wrapped_sum(vec, 4);
    double sum_kah = kahan_sum(vec, 4);
    double sum_kbn = kbn_sum(vec, 4);

    printf("Subtask A -----------------------\n");
    printf(" for loop sum:    %.17g\n", sum_for);
    printf(" GSL sum:         %.17g\n", sum_gsl);
    printf(" Kahan sum:       %.17g\n", sum_kah);
    printf(" (extra) KBN sum: %.17g\n", sum_kbn);


    // subtask b -----------------------
    srand( time(NULL) );
    int n = 10000000;
    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));

    fill_gaussian(x, n);
    fill_gaussian(y, n);

    // DAXPY (inplace operation on y)
    daxpy(1.0, x, y, n);

    // IDEA: compute mean and std of resulting y,
    // which should be Gaussian with mean=0 and std=sqrt(2)

    double res_mean = compute_mean(y, n);
    double res_stdv = compute_std(y, n, res_mean);
    printf("\nSubtask B -----------------------\n");
    printf(" mean = %.5f, std = %.5f\n", res_mean, res_stdv);

    double expected_mean = 0.0;
    double expected_stdv = sqrt(2.0);
    assert(fabs(res_mean - expected_mean) < 1e-2);
    assert(fabs(res_stdv - expected_stdv) < 1e-2);
    printf(" checks passed!\n");

    free(x);
    free(y);

    return 0;
}
