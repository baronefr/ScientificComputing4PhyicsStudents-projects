#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./fileio.h"

#define STRLEN 512

void generate_vector(int N, double **v, double init_value) {
    // allocate memory for vector, if not yet allocated
    if (*v == NULL) {
        printf("Allocating memory for vector.\n");
        *v = (double*)malloc(N*sizeof(double));
        if (*v == NULL) {
            printf("Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
    }
    // else: we assume the memory is already allocated and
    // that it fits N double values

    // initialize vector
    for (int i=0; i<N; i++) (*v)[i] = init_value;
}

int main(int argc, char** argv) {

    // parse command line arguments
    if (argc != 3) {
        printf("Wrong input pattern: %s <N> <fname_prefix>\n", argv[0]);
        return EXIT_FAILURE;
    }
    int N = atoi(argv[1]);
    const char *fname_prefix = argv[2];
    printf("Generating vectors with N=%d elements on files %s*\n", N, fname_prefix);

    // variables
    double *pvect = NULL;
    char filename[STRLEN];
    
    // first handle X...
    generate_vector(N, &pvect, 0.1);

    snprintf(filename, sizeof(filename), "%s_N%d_x.dat", fname_prefix, N);
    dump_vector_binary(N, filename, pvect);

    // ... then Y
    generate_vector(N, &pvect, 7.1);

    snprintf(filename, sizeof(filename), "%s_N%d_y.dat", fname_prefix, N);
    dump_vector_binary(N, filename, pvect);

    // closing operations
    free(pvect);

    return EXIT_SUCCESS;
}
