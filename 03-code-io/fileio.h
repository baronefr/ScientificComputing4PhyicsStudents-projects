#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

void read_vector_binary(int N, const char *fname, double **vector) {
    FILE *file = fopen(fname, "rb");
    if (file == NULL) {
        printf("Error: cannot open file <%s>\n", fname);
        exit(EXIT_FAILURE);
    }

    *vector = (double*)malloc(N*sizeof(double));
    if (*vector == NULL) {
        printf("Memory allocation failed\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fread(*vector, sizeof(double), N, file);
    fclose(file);
}

void dump_vector_binary(int N, const char *fname, double *vect) {
    FILE *file;

    printf("Writing file <%s> ...", fname);
    file = fopen(fname, "wb");
    if (file == NULL) {
        printf(" error: cannot open file.\n");
        free(vect);
        exit(EXIT_FAILURE);
    } else fwrite(vect, sizeof(double), N, file);
    fclose(file);
    printf(" done.\n");
}
