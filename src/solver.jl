export solve, solve!


function solve(model::CNLSModel; silent::Bool=false, max_iter::Int64 = 100, scaling::Bool=false)
    ε = eps(eltype(model.starting_point))
    sqr_ε = sqrt(ε)
    exit_code, x_opt, f_opt = enlsip(model.starting_point, model.residuals, model.constraints, model.nb_parameters, model.nb_residuals, model.nb_eqcons, model.nb_cons, 
    verbose=!silent, scaling=scaling, MAX_ITER=max_iter, ε_rel = sqr_ε, ε_x = sqr_ε, ε_c = sqr_ε)
    sol = CnlsResult(exit_code >0, x_opt, f_opt)
    return sol
end

"""
    solve!(model)

Once a [`CnlsModel`](@ref) has been instantiated, this function solves the optimzation problem associated by using the method implemented in Enlsip.


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
function solve!(model::CnlsModel; silent::Bool=false, max_iter::Int64=100, scaling::Bool=false)

    # Relative precision
    ε = eps(eltype(model.starting_point))
    sqr_ε = sqrt(ε)
    
    # Evaluation functions
    residuals_evalfunc = (model.jacobian_residuals === nothing ? ResidualsFunction(model.residuals) : ResidualsFunction(model.residuals, model.jacobian_residuals))

    if all(!isfinite,vcat(model.x_low, model.x_upp))
        constraints_evalfunc = instantiate_constraints_wo_bounds(model.eq_constraints, model.jacobian_eqcons, model.ineq_constraints, model.jacobian_ineqcons)
    else
        constraints_evalfunc = instantiate_constraints_w_bounds(model.eq_constraints, model.jacobian_eqcons, model.ineq_constraints, model.jacobian_ineqcons, model.x_low, model.x_upp)
    end

    nb_constraints = total_nb_constraints(model)

    # Call ENLSIP-Julia
    exit_code, x_opt, f_opt = enlsip(model.starting_point, residuals_evalfunc, constraints_evalfunc, model.nb_parameters, model.nb_residuals, model.nb_eqcons, nb_constraints,
    verbose=!silent, scaling=scaling, MAX_ITER=max_iter, ε_rel = sqr_ε, ε_x = sqr_ε, ε_c = sqr_ε)

    # Update of solution related fields in model
    model.status_code = convert_exit_code(exit_code)
    model.sol = x_opt
    model.obj_value = f_opt
    return
end