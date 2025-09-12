
import numpy as np


def stats_errors(A, Arec, eps=1e-200):
    """
    Statistics of reconstruction errors between original matrix A and its reconstruction Arec.
    """
    diff = Arec - A
    abs2 = diff.ravel() ** 2 # squared absolute errors
    rel = diff / (np.abs(A) + eps)
    rel2 = rel.ravel() ** 2 # squared relative errors
    return {
        "rmse_abs": np.sqrt( np.real(abs2.mean()) ),
        "rmed_abs": np.sqrt( np.real(np.median(abs2)) ),
        "rmse_rel": np.sqrt( np.real(rel2.mean()) ),
        "rmed_rel": np.sqrt( np.real(np.median(rel2)) ),
    }


def cFFT_from_r2c(R, n):
    """
    Reconstructs cFFT from rFFT output.
    F[i,j] = conj(F[-i mod n, -j mod n]).
    """
    m = n // 2 + 1
    F = np.zeros((n, n), dtype=np.complex128)

    is_even_columns = (n % 2 == 0)

    # copy half of matrix
    F[:, :m] = R
    # fill remaining columns using F[i,j] = conj(F[-i mod n, -j mod n])
    for k in range(n):
        k_sym = (-k) % n
        for l in range(1, m):
            if is_even_columns and (l == n // 2):
                # small optimization: skip Nyquist column for even n (self-symmetric)
                continue
            l_sym = n - l
            F[k, l_sym] = np.conj(R[k_sym, l])
    return F


if __name__ == "__main__":
    n = 1000
    A = np.random.normal(loc=1,scale=1,size=(n,n))
    A_complex = A.astype(np.complex128)

    # Questions 1 & 2: c2c FFT and inverse ---
    C = np.fft.fft2(A_complex)
    A_rec_c2c = np.fft.ifft2(C).real
    stats_c2c = stats_errors(A_complex, A_rec_c2c)

    print("Questions 1&2 (c2c):")
    print(f"| RMSE abs  : {stats_c2c['rmse_abs']:.3e}")
    print(f"| RMED abs  : {stats_c2c['rmed_abs']:.3e}")
    print(f"| RMSE rel  : {stats_c2c['rmse_rel']:.3e}")
    print(f"| RMED rel  : {stats_c2c['rmed_rel']:.3e}\n")

    # Questions 3 & 4: r2c FFT and inverse ---
    R = np.fft.rfft2(A)
    A_rec_r2c = np.fft.irfft2(R, s=A.shape)
    stats_r2c = stats_errors(A, A_rec_r2c)

    print("Questions 3&4 (r2c):")
    print(f"|  RMSE abs  : {stats_r2c['rmse_abs']:.3e}")
    print(f"|  RMED abs  : {stats_r2c['rmed_abs']:.3e}")
    print(f"|  RMSE rel  : {stats_r2c['rmse_rel']:.3e}")
    print(f"|  RMED rel  : {stats_r2c['rmed_rel']:.3e}\n")

    # Question 5: machine precision ---
    eps = np.finfo(float).eps
    print("Question 5:")
    print(f"| Machine precision eps: {eps}\n")

    # Question 6: C[0,0] and R[0,0] ---
    C00 = C[0, 0]
    R00 = R[0, 0]

    print("Question 6:")
    print(f"| C[0,0]: {C00}")
    print(f"| R[0,0]: {R00}")
    print(">  Note: C[0,0] = R[0,0] = sum of all elements in A.")
    print(f">    Indeed, sum(A) = {np.sum(A)}\n")

    # Question 7: 6x6 test (bonus) ---
    n_small = 6
    A6 = np.random.normal(loc=1, scale=1, size=(n_small, n_small))
    A6_complex = A6.astype(np.complex128)

    C6 = np.fft.fft2(A6_complex)
    R6 = np.fft.rfft2(A6)
    C6_from_R6 = cFFT_from_r2c(R6, n_small)

    print("Question 7 (6x6):")
    assert np.allclose(C6, C6_from_R6), "error: C6 and C6_from_R6 are not close"  # should be True
    print(f"|  test passed: C6 and C6_from_R6 are close")
