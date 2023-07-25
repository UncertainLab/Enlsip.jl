export solve

"""
    solution = solve(model)

Once a [`EnlsipModel`](@ref) has been instantiated, this function solves the optimzation problem associated by using the method implemented in Enlsip.

This function returns an object of type `EnlsipSolution`.
"""
function solve(model::EnlsipModel; silent::Bool=false, max_iter::Int64 = 100, scaling::Bool=false)
    ε = eps(eltype(model.starting_point))
    sqr_ε = sqrt(ε)
    sol = enlsip(model.starting_point, model.residuals, model.constraints, model.nb_parameters, model.nb_residuals, model.nb_eqcons, model.nb_cons, 
    verbose=!silent, scaling=scaling, MAX_ITER=max_iter, ε_rel = sqr_ε, ε_x = sqr_ε, ε_c = sqr_ε)
    return sol
end