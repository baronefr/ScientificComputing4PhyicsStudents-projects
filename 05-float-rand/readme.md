# task 05 handout

## Subtask A

All the three proposed algorithms to sum the elements of the array fail.

The trivial for-loop sum fails because:
- `1.0 + 1.0e16 = 1.0e16`: the 1.0 is lost because of double precision rounding
- `1.0e16 - 1.0e16 = 0.0`, correctly
- `0.0 - 0.5 = -0.5`, correctly from the numerical point of view, but returns the wrong algebraic result (we expect 0.5)

The GSL routine is expected to fail, as it implements the same for-loop under the hood.

Kahan fails as well. The reason is that the comulative variable `c` keeps track of the **tiny bits** lost when adding small and large numbers, but this is not effective for this case.
> In the second iteration of Kahan's sum, the loop at `i=1` starts with `y=1, t=1, c=0, sum=1`.
> The instruction `t = sum + y` is affected by catastrofic cancelation, since `sum=1, y=1e16`. The following instruction, `c = (t - sum) - y` is affected by catastrofic cancelation too, since it expands to `(1e16-1) - 1e16`.

### extra

I implement [Neumaier's variant](https://en.wikipedia.org/wiki/Kahan_summation_algorithm) of the Kahan algorithm, which is able to fix the sum.

## Subtask B

To check for correctness, I test the statistical properties of the resulting sum.
