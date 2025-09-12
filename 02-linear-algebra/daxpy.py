
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
