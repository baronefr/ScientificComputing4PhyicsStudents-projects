#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

void read_vector_hdf5(int N, const char *fname, const char *dataset_name, double **vector) {
    hid_t file_id, dataset_id;
    herr_t status;

    // IO
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);

    // allocate memory
    *vector = (double*)malloc(N*sizeof(double));
    if (*vector == NULL) {
        printf("Memory allocation failed\n");
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        exit(EXIT_FAILURE);
    }

    // read from file
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *vector);

    // close
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    if (status < 0) {
        printf("Error reading dataset %s from file %s\n", dataset_name, fname);
        free(*vector);
        exit(EXIT_FAILURE);
    }
}

void dump_vector_hdf5(int N, const char *fname, const char *dataset_name, double *vect) {
    hid_t file_id, dataset_id, dataspace_id;
    hsize_t dims[1] = {N};

    printf("Writing file <%s> ... ", fname);

    // create a new file & dataset
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // write
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vect);

    // closing operations
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);

    printf("done.\n");
}
