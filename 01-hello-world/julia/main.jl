using ArgParse
using Printf
using DelimitedFiles

const DEF_OUTPUT_DIR = "output"

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--help, -h" => help="produce help message"
        "--outdir, -o" => ArgParse.StoreValue{String}; default=DEF_OUTPUT_DIR; help="output directory"
    end

    parsed_args = parse_args(s)
    
    if parsed_args["help"]
        println(s)
        return
    end

    OUTPUT_DIR = parsed_args["outdir"]
    OUTPUT_FILE = joinpath(OUTPUT_DIR, "output_julia.txt")

    # Create the output directory
    mkpath(OUTPUT_DIR)

    # Define the input values
    a = 3.0
    n = 20
    x = fill(1.0, n)
    y = fill(4.0, n)

    # Calculate z = a*x + y
    z = a .* x .+ y

    # Write the result to a text file
    open(OUTPUT_FILE, "w") do io
        for val in z
            println(io, val)
        end
    end
end

main()