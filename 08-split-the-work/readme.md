# task 08 handout

No comments or questions to answer, but I would have implemented the chunked daxpy differently. For instance, the elements that not fit into a chunck can be handled outside the chunk loop, to avoid the if statement at every iteration.
