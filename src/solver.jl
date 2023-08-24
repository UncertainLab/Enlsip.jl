export solve

"""
    solution = solve(model)

Once a [`CNLSModel`](@ref) has been instantiated, this function solves the optimzation problem associated by using the method implemented in Enlsip.

This function returns an object of type [`CNLSResult`](@ref).

The following optionnal arguments can be provided:

* `silent::Bool` 
    
    - Set to `true` if one wants the algorithm to print details about the iterations and termination of the solver

    - Default value is set to `false`

* `max_iter::Int` 

    - Maximum number of iterations allowed

    - Default value is set to `100`

* `scaling::Bool` 

    - Set to `true` if one wants the algorithm to work with a constraints jacobian matrix whose rows are scaled (i.e. all constraints gradients vectors are scaled)

    - Default value is set to `false`
"""
function solve(model::CNLSModel; silent::Bool=false, max_iter::Int64 = 100, scaling::Bool=false)
    ε = eps(eltype(model.starting_point))
    sqr_ε = sqrt(ε)
    sol = enlsip(model.starting_point, model.residuals, model.constraints, model.nb_parameters, model.nb_residuals, model.nb_eqcons, model.nb_cons, 
    verbose=!silent, scaling=scaling, MAX_ITER=max_iter, ε_rel = sqr_ε, ε_x = sqr_ε, ε_c = sqr_ε)
    return sol
end