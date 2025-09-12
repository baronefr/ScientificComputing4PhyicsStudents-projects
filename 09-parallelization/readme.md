# task 09 handout

No comments or questions to answer.
I slightly don't like the need to broadcast the vector `d` to all processes (see line 171):
```C
MPI_Bcast(d_ref, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
```
The best solution would be to broadcast to each MPI process only the chunk to sum.
