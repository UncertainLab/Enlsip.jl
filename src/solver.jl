export solve!, print_cnls_model

"""
    solve!(model)

Once a [`CnlsModel`](@ref) has been instantiated, this function solves the optimzation problem associated by using the method implemented in the `Enlsip` solver.


Keywords arguments:

* `silent::Bool` 
    
    - Set to `false` if one wants the algorithm to print details about the iterations and termination of the solver

    - Default value is `true`, i.e. by default, there is no output. If one wants to print those information afert solving, the [`print_cnls_model`](@ref) method 
    can be called.

* `max_iter::Int` 

    - Maximum number of iterations allowed

    - Default value is set to `100`

* `scaling::Bool` 

    - Set to `true` if one wants the algorithm to work with a constraints jacobian matrix whose rows are scaled (i.e. all constraints gradients vectors are scaled)

    - Default value is set to `false`
"""
function solve!(model::CnlsModel; silent::Bool=true, max_iter::Int=100, scaling::Bool=false)

    # Relative precision
    ε = eps(eltype(model.starting_point))
    sqr_ε = sqrt(ε)
    
    # Internal scaling
    model.constraints_scaling = scaling
    
    # Evaluation functions
    residuals_evalfunc = (model.jacobian_residuals === nothing ? ResidualsFunction(model.residuals) : ResidualsFunction(model.residuals, model.jacobian_residuals))

    if all(!isfinite,vcat(model.x_low, model.x_upp))
        constraints_evalfunc = instantiate_constraints_wo_bounds(model.eq_constraints, model.jacobian_eqcons, model.ineq_constraints, model.jacobian_ineqcons)
    else
        constraints_evalfunc = instantiate_constraints_w_bounds(model.eq_constraints, model.jacobian_eqcons, model.ineq_constraints, model.jacobian_ineqcons, model.x_low, model.x_upp)
    end

    nb_constraints = total_nb_constraints(model)

    # Call the ENLSIP solver
    exit_code, x_opt, f_opt, enlsip_info = enlsip(model.starting_point, residuals_evalfunc, constraints_evalfunc, model.nb_parameters, model.nb_residuals, model.nb_eqcons, nb_constraints,
    scaling=scaling, MAX_ITER=max_iter, ε_rel = sqr_ε, ε_x = sqr_ε, ε_c = sqr_ε)

    # Update of solution related fields in model
    model.model_info = enlsip_info
    model.status_code = convert_exit_code(exit_code)
    model.sol = x_opt
    model.obj_value = f_opt

    !silent && print_diagnosis(model)
    return
end

"""
    print_cnls_model(model,io)

One can call this function to print information about an instance `model` (see [`CnlsModel`](@ref)). 

If `model` has just been instantiated but not solved, it will print general information about the model, such as the dimensions of the residuals, parameters and constraints.

After calling the [`solve!`](@ref) method, the output will be enhanced with details about the iterations performed during the execution of the algorithm.

The following info are also printed:

* number of iterations

* total number of function and Jacobian matrix evaluations for both residuals and contraints

* solving time in seconds

* value of the objective function found by the algorithm

* termination status (see [`status`](@ref))
"""
function print_cnls_model(model::CnlsModel, io::IO=stdout)
    model_status = status(model)
    if model_status == :unsolved
        print_initialized_model(model, io)
    else
        print_diagnosis(model, io)
    end
end