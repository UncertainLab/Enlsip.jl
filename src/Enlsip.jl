module Enlsip

# Packages

using LinearAlgebra, Polynomials
using Formatting, Printf

# include source files
for f in ["structures", "enlsip_functions", "solver"]
    include("./$f.jl")
end

end # module