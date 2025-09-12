# task 06 handout

The text in the boxes represents the outputs from `main.py`.

## 1+2: c2c FFT

This question asks us to perform c2c FFT of matrix $A$ and reconstruct the original matrix by inverting the FFT ($A_r$). The Python script then computes statistical indicators of the reconstruction accuracy:
$$RMSE_{abs} = \sqrt{\mathrm{mean}((A-A_r)^2)}$$
$$RMED_{abs} = \sqrt{\mathrm{median}((A-A_r)^2)}$$
$$RMSE_{rel} = \sqrt{\mathrm{mean}(((A-A_r)/A)^2)}$$
$$RMED_{rel} = \sqrt{\mathrm{median}(((A-A_r)/A)^2)}$$

```bash
Questions 1&2 (c2c):
| RMSE abs  : 5.360e-16
| RMED abs  : 4.198e-16
| RMSE rel  : 1.687e-13
| RMED rel  : 3.570e-16
```

The FFT is accurate within machine precision ($\sim 2.22*10^{-16}$).

## 3+4: r2c FFT

```bash
Questions 3&4 (r2c):
|  RMSE abs  : 6.102e-16
|  RMED abs  : 4.441e-16
|  RMSE rel  : 1.487e-13
|  RMED rel  : 3.970e-16
```

## 5: machine precision

I answered in the previous sections. FFT from scipy is accurate.

## 6: first element of FFT

The element $[0,0]$ of FFT is equal to the sum of all elements in the matrix. I have validated this by explicit computation.

```bash
Question 6:
| C[0,0]: (1001472.3866234156+0j)
| R[0,0]: (1001472.3866234156+0j)
>  Note: C[0,0] = R[0,0] = sum of all elements in A.
>    Indeed, sum(A) = 1001472.3866234155
```

## 7 (bonus): cFFT from rFFT

The FFT of a real matrix is hermitian-symmetric: the negative frequency terms are just the complex conjugates of the corresponding positive-frequency terms, and the negative-frequency terms are therefore redundant. 
The rFFT routine from scipy returns only the positive-frequency terms in form of a triangular matrix.
Knowing this, the output of a complex FFT (which instead computes explicitly all the frequency terms) can be obtained from the complex conjugate of the rFFT matrix:
$$F[i,j] = conj(F[-i mod n, -j mod n])\quad .$$
