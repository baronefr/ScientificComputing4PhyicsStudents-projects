# task 02 handout

submitted by Barone Francesco Pio

## Q&A

1) Did you find any problems in running the codes for some $N$. If so, do you have an idea why?

    > No, I had no problem. But I might guess that if you had one, it would be to allocate the memory for some size-$N$ matrices. $10000^2 \cdot 64 \mathrm{bits} \simeq 744 \mathrm{MB}$ of memory. For $N=100000$ it would be already $\sim 80\mathrm{GB}$.

2) Were you able to test correctly the sum and product of points 1-3? If so, how? If not, what was the problem?

    > I was. The catch here is to check up to some reasonable precision, because a naive check like `d[i] == 7.4` would fail due to machine precision issues. Also, the question of point 3 is obviously wrong, since for $c_{ij}=a_{ik}b_{kj}$, the result would be $c_{ij}=0.3\cdot7.1\cdot N$ - which $=21.3$ only for $N=10$.

## Source codes

I append in the following the source code for my tests, in C and python.

### axpy

```C
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// example constants
#define EXAMPLE_ARRAY_SIZE 100000000 // also test 10, 1000000 
#define EXAMPLE_EXPECTED_RESULT 7.4
#define DOUBLE_EPSILON 1e-14


void daxpy_not_inplace(int n, double a, double *px, double *py, double *pd) {
    // computes pd = a*px + py
    // where int n is the size of arrays px, py, pd
    int i, m;

    // return on trivial cases
    if (n<=0) return;
    if (a==0.0) {
        for (i=0; i<n; i++) pd[i] = 0.0;
    }

    // clean-up loop
    m = n%4;
    if (m!=0) {
        for (i=0; i<m; i++) {
            pd[i] = py[i] + a * px[i];
        }
        if (n<4) return;
    }
    // loop unrolling in units of 4
    // (assuming 256-bit registers)
    for (i=m; i<n; i+=4) {
        pd[i]   = py[i] + a*px[i];
        pd[i+1] = py[i+1] + a*px[i+1];
        pd[i+2] = py[i+2] + a*px[i+2];
        pd[i+3] = py[i+3] + a*px[i+3];
    }

}


int main() {
    double a;
    double *px, *py, *pd;
    int i = 0;
    bool test = true;

    // allocating arrays
    printf("allocating...\n");
    px = (double*)malloc(EXAMPLE_ARRAY_SIZE * sizeof(double));
    py = (double*)malloc(EXAMPLE_ARRAY_SIZE * sizeof(double));
    pd = (double*)malloc(EXAMPLE_ARRAY_SIZE * sizeof(double));

    if (px == NULL || py == NULL || pd == NULL) {
        printf("[error] memory allocation failed\n");
        return 1;
    }

    // initialize arrays & constants for this test
    for (i = 0; i < EXAMPLE_ARRAY_SIZE; i++) {
        px[i] = 0.1;
        py[i] = 7.1;
    }
    a = 3;

    // daxpy (note: output on pd)
    printf("computing...\n");
    daxpy_not_inplace(EXAMPLE_ARRAY_SIZE, a, px, py, pd);

    // testing the result
    printf("checking result...");
    for (i = 0; i < EXAMPLE_ARRAY_SIZE; i++)
        if (fabs(pd[i] - EXAMPLE_EXPECTED_RESULT) >= DOUBLE_EPSILON) break;

    // printing the result
    if (i < EXAMPLE_ARRAY_SIZE) {
        printf(" test failed!\n");
        printf("  expected result is %f, but for instance d[%d]=%f\n", EXAMPLE_EXPECTED_RESULT, i, pd[i]);
    }
    else printf(" all elements are correct\n");

    // free allocated memory
    free(px);
    free(py);
    free(pd);

    return 0;
}
```

```python
import numpy as np

N = 100000000

a = 3
x = 0.1*np.ones(N)
y = 7.1*np.ones(N)

d = a*x + y

# check
test = np.allclose(d, 7.4)
print(f"test result is: {test}")
assert test
```

### matmul

```C
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void matmul(double *A, double *B, double *C, int n) {
    // naive way to do matmul ...
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            C[i*n + j] = 0.0;
            for (int k = 0; k < n; k++)
                C[i*n + j] += A[i*n + k] * B[k*n + j];
        }
}

bool allclose(double *C, double expected, int n, double epsilon) {
    for (int i = 0; i < n*n; i++)
        if (fabs(C[i] - expected) >= epsilon) return false;
    return true;
}

int main() {
    int N[] = {10, 100, 10000};
    int num_tests = sizeof(N) / sizeof(N[0]);
    double *A, *B, *C;
    int nn;
    bool test;

    for (int t = 0; t < num_tests; t++) {
        nn = N[t];
        printf("testing for N=%d ...\n", nn);

        A = (double*)malloc(nn * nn * sizeof(double));
        B = (double*)malloc(nn * nn * sizeof(double));
        C = (double*)malloc(nn * nn * sizeof(double));

        if (A == NULL || B == NULL || C == NULL) {
            printf("[error] memory allocation failed\n");
            return 1;
        }

        for (int i = 0; i < nn * nn; i++) {
            A[i] = 0.3;
            B[i] = 7.1;
        }

        matmul(A, B, C, nn);

        printf(" checking the result... ");
        test = allclose(C, 0.3 * 7.1 * nn, nn, 1e-8);
        if (test) printf(" passed\n");
        else {
            printf(" failed\n");
            printf("  C[0] = %.15f but expected %.15f\n", C[0], 0.3 * 7.1 * nn);
            return 1;
        }

        free(A);
        free(B);
        free(C);
    }

    return 0;
}
```

```python

import numpy as np

N = [ 10, 100, 10000 ]

for nn in N:
    print(f"testing for N={nn}...")
    
    A = np.ones((nn,nn))*0.3
    B = np.ones((nn,nn))*7.1

    C = A @ B

    # check
    test = np.allclose(C, 0.3*7.1*nn)
    print(f" test result for N={nn} is {test}")
    assert test
```
