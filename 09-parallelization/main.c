#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <assert.h>


// quick wall time
double wall_time() {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

// fill with random doubles in [-1,1]
void fill_random(double *v, int n) {
    for (int i = 0; i < n; i++) {
        v[i] = 2.0 * rand() / (double)RAND_MAX - 1.0;
    }
}



// DAXPY -------------

void daxpy(double *d, double a, double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        d[i] = a * x[i] + y[i];
    }
}

void daxpy_openmp(double *d, double a, double *x, double *y, int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) d[i] = a * x[i] + y[i];
}

void daxpy_mpi(double *d, double a, double *x, double *y, int n, int rank, int nprocs, MPI_Comm comm) {
    int chunk = (n + nprocs - 1) / nprocs;
    int start = rank * chunk;
    int end   = (start + chunk > n) ? n : (start + chunk);

    for (int i = start; i < end; i++) {
        d[i] = a * x[i] + y[i];
    }

    // gather results on root (rank 0)
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, d, chunk, MPI_DOUBLE, comm);
}


// reduction -------------

double reduce_openmp(double *d, int n) {
    double s = 0.0;
    #pragma omp parallel for reduction(+:s)
    for (int i = 0; i < n; i++) s += d[i];
    return s;
}

double reduce_mpi(double *d, int n, MPI_Comm comm) {
    // get rank and nprocs (again, could be passed as arguments)
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    double local = 0.0, global = 0.0;

    int chunk = (n + nprocs - 1) / nprocs;
    int start = rank * chunk;
    int end   = (start + chunk > n) ? n : (start + chunk); // limit upper limit to n

    for (int i = start; i < end; i++) 
        local += d[i];

    MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    return global;
}



// main -------------

int main(int argc, char **argv) {
    int N = 1000;  // default size, small for debugging
    if (argc > 1) N = atoi(argv[1]);

    srand(12345);

    // initialize MPI
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // create data for daxpy (random)
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *d = malloc(N * sizeof(double));
    double *d_ref = malloc(N * sizeof(double));
    double a = 2.0;
    double sum_serial = 0.0;
    fill_random(x, N);
    fill_random(y, N);

    // 0: serial daxpy, for reference
    if (rank == 0) {
        double tic = wall_time();
        daxpy(d_ref, a, x, y, N);
        double toc = wall_time();
        printf("[serial] time = %.6f s\n", toc - tic);
        for (int i = 0; i < N; i++) 
            sum_serial += d_ref[i];
    }

    // 1: OpenMP daxpy
    if (rank == 0) {
        double tic = wall_time();
        daxpy_openmp(d, a, x, y, N);
        double toc = wall_time();
        printf("[OpenMP] time = %.6f s\n", toc - tic);

        // check
        for (int i = 0; i < N; i++) {
            if (fabs(d[i] - d_ref[i]) > 1e-12) {
                printf("Mismatch at %d!\n", i);
                break;
            }
        }
        printf(" [test] daxpy openmp vs serial passed\n");
    }

    // let's reset d to random values
    //fill_random(d, N);


    // 2: MPI daxpy
    MPI_Barrier(MPI_COMM_WORLD);
    double tic = wall_time();
    daxpy_mpi(d, a, x, y, N, rank, nprocs, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    double toc = wall_time();
    if (rank == 0) {
        printf("[MPI]    time = %.6f s\n", toc - tic);

        // check
        for (int i = 0; i < N; i++) {
            if (fabs(d[i] - d_ref[i]) > 1e-12) {
                printf("MPI mismatch at %d!\n", i);
                break;
            }
        }
        printf(" [test] daxpy mpi vs serial passed\n");
    }



    // (optional) reduction
    if (rank == 0) {
        printf("\n(optional) reduction -----\n");
        double sum_omp = reduce_openmp(d_ref, N);
        printf("[serial] sum = %.10f\n", sum_serial);
        printf("[OpenMP] sum = %.10f\n", sum_omp);
        assert(fabs(sum_serial - sum_omp) < 1e-10);
        printf(" [test] reduce openmp vs serial passed\n");
    }

    // FIXME: broadcast d_ref
    MPI_Bcast(d_ref, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double sum_mpi = reduce_mpi(d_ref, N, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("[MPI]    sum = %.10f\n", sum_mpi);
        assert(fabs(sum_serial - sum_mpi) < 1e-10);
        printf(" [test] reduce mpi vs serial passed\n");
    }

    // cleaning
    free(x);
    free(y);
    free(d);
    free(d_ref);
    MPI_Finalize();

    return 0;
}
