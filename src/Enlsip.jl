module Enlsip

# Packages

using LinearAlgebra, Polynomials, ForwardDiff, Printf

# include source files
for f in ["structures","cnls_model", "enlsip_functions", "solver"]
    include("./$f.jl")
end

end # module