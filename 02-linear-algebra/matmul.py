
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
