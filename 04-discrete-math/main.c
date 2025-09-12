#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

// abs ABS_TOLERANCE for the check of the integral convergence (in the known case)
#ifndef ABS_TOLERANCE
    #define ABS_TOLERANCE 1e-8
#endif

// if defined, the sampled points will be saved to a file
//#define FILE_OUTPUT "samples.txt"

// f(x)
double f(double x) {
    return exp(x)*cos(x);
}

double f_primitive(double x) {
    return (exp(x)*(cos(x) + sin(x)))/2.0;
}

// integrate using the trapezoidal rule
// see: https://en.wikipedia.org/wiki/Trapezoidal_rule
double integral_by_trapezoid(int N, double x_inf, double x_sup) {
    // memory-conservative approach: do not store the sampled points, but loop in the sum
    double integral;
    double x;

    if(N < 3) {
        printf("N must be greater than 2\n");
        exit(EXIT_FAILURE);
    }

    integral = (f(x_inf) + f(x_sup))/2.0;

    for (int i = 1; i < N - 1; i++) {
        x = x_inf + ( (i*(x_sup-x_inf))/(N) );
        integral += f(x);
    }

    integral = (integral*(x_sup - x_inf)) / N;

    return integral;
}

// evaluate if two doubles are close
bool isclose(double C, double expected, double epsilon) {
    if (fabs(C - expected) >= epsilon) return false;
    return true;
}


int main(int argc, char** argv) {

    // parse command line arguments
    if (argc != 4) {
        printf("Wrong input pattern: <N> <low_limit> <up_limit>\n");
        return EXIT_FAILURE;
    }

    int N;
    double integral, x_inf, x_sup, expected, rel_error;

    N = atoi(argv[1]); // NOTE: N is number of trapezoids, so N+1 points will be sampled
    x_inf = atof(argv[2]);
    x_sup = atof(argv[3]);
    printf("Integration parameters: N=%d, x_inf=%f, x_sup=%f\n", N, x_inf, x_sup);

    #ifdef FILE_OUTPUT
        double x, fx;
        printf("saving samples to file ... ");
        
        FILE *file = fopen(FILE_OUTPUT, "w");

        if (file == NULL) {
            perror("error opening file\n");
            return EXIT_FAILURE;
        }
        
        // NOTE: sampling N+1 points
        for (int i = 0; i <= N; i++) {
            x = x_inf + (i*(x_sup-x_inf))/N;
            fx = f(x);
            fprintf(file, "%.16f\t%.16f\n", x, fx);
        }
        fclose(file);
        printf("done\n");
    #endif

    // compute int of f(x)
    printf("computing ... ");
    integral = integral_by_trapezoid(N, x_inf, x_sup);
    printf("done\n");

    printf("\n> integration interval: (%.8f,%.8f)\n> int f dx = %.16f\n", x_inf, x_sup, integral);

    // check the result
    if ( isclose(x_inf, 0.0, 1e-6) & isclose(x_sup, M_PI/2, 1e-6) ) {
        // this is the known case, requested in the exercise
        expected = (exp(M_PI/2)-1.0)/2.0;
    } else {
        // the primitive of f(x) is known ... so we can check in general cases too
        expected = f_primitive(x_sup) - f_primitive(x_inf);
    }

    rel_error = (integral/expected) - 1;
    printf("  relative error = %.10e\n", rel_error);

    if ( isclose(integral, expected, ABS_TOLERANCE) ) {
        printf("Test passed up to abs precision %.3e\n", ABS_TOLERANCE);
        return EXIT_SUCCESS;
    } else {
        printf("Test failed up to abs precision %.3e\n", ABS_TOLERANCE);
        printf("| expected: %.16f, but integration returns %.16f\n", expected, integral);
        
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
