#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./parser.h"
#include "./fileio.h"

#define STRLEN 512

void daxpy(int n, double a, double *px, double *py) {
    // computes py = a*px + py
    // where int n is the size of arrays px, py
    int i, m;

    // return on trivial cases
    if (n<=0) return;
    if (a==0.0) {
        for (i=0; i<n; i++) py[i] = 0.0;
    }

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
    int N = 0;
    double a = 0.0;
    double *px, *py;
    char fname_x[STRLEN], fname_y[STRLEN], fname_prefix[STRLEN], filename[STRLEN];

    // parsing inputs from conf file
    config_option_t co;
    if ((co = read_config_file("daxpy.conf")) == NULL) {
        printf("read_config_file() error");
        return EXIT_FAILURE;
    }
    printf("parsed from .conf:\n");
    while(1) {
        printf(" | %s -> %s\n", co->key, co->value);

        if (strcmp(co->key, "fname_x") == 0) {
            strcpy(fname_x, co->value);
        } else if (strcmp(co->key, "fname_y") == 0) {
            strcpy(fname_y, co->value);
        } else if (strcmp(co->key, "N") == 0) {
            N = atoi(co->value);
        } else if (strcmp(co->key, "a") == 0) {
            a = atof(co->value);
        } else if (strcmp(co->key, "prefix") == 0) {
            strcpy(fname_prefix, co->value);
        }
        
        if (co->prev != NULL) co = co->prev;
        else break;
    }
    printf("N=%d x=\"%s\" y=\"%s\" a=%f prefix=\"%s\"\n", N, fname_x, fname_y, a, fname_prefix);

    // reading x and y
    read_vector_binary(N, fname_x, &px);
    read_vector_binary(N, fname_y, &py);

    daxpy(N, a, px, py); // result is in py
    printf("[info] look at first result: d[0] = %f\n", py[0]);

    // dump to file
    snprintf(filename, sizeof(filename), "%s_N%d_d.dat", fname_prefix, N);
    dump_vector_binary(N, filename, py);
    free(px);
    free(py);

    return EXIT_SUCCESS;
}
