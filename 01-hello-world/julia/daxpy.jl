
# explicitly unrolled loop daxpy 
function daxpy(n::Int, a::Float64, px::Vector{Float64}, py::Vector{Float64})
    # return on trivial cases
    if n <= 0 || a == 0.0
        return
    end

    # clean-up loop
    m = n%4
    for i in 1:m
        py[i] += a*px[i]
    end
    # loop unrolling in units of 4
    # (assuming 256-bit registers)
    for i in m+1:4:n-3
        py[i]   += a*px[i]
        py[i+1] += a*px[i+1]
        py[i+2] += a*px[i+2]
        py[i+3] += a*px[i+3]
    end

end

# not unrolled loop daxpy 
function daxpy_nounrolling(n::Int, a::Float64, px::Vector{Float64}, py::Vector{Float64})
    # return on trivial cases
    if n <= 0 || a == 0.0
        return
    end

    # no loop unrolling here
    for i in 1:n
        py[i] += a*px[i]
    end

end


# Allocate and initialize vectors
EXAMPLE_ARRAY_SIZE = 1000000
EXAMPLE_A = 3.14
px = rand(EXAMPLE_ARRAY_SIZE)
py = rand(EXAMPLE_ARRAY_SIZE)

# daxpy (note: in-place on py)
daxpy(EXAMPLE_ARRAY_SIZE, EXAMPLE_A, px, py)

# Print the first few results
println("Printing first 10 elements:")
println( py[1:10] )


using BenchmarkTools

println("benchmarking unrolled routine:")
display( @benchmark daxpy($EXAMPLE_ARRAY_SIZE, $EXAMPLE_A, $px, $py) )

println("benchmarking not-unrolled routine:")
display( @benchmark daxpy_nounrolling($EXAMPLE_ARRAY_SIZE, $EXAMPLE_A, $px, $py) )

# comments ----
#
# I expect the not-unrolled routine to work as efficient as the explicitly unrolled one,
# maybe even better. The reason is that Julia compiler is optimizing it nevertheless.
# Check `@code_llvm debuginfo=:none` to be sure about this.
