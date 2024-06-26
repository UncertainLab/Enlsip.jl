export solve!, print_cnls_model

"""
    solve!(model{T})

Once a [`CnlsModel`](@ref) has been instantiated, this function solves the optimzation problem associated by using the method implemented in the `Enlsip` solver. 

Options:

* `silent::Bool` 
    
    - Set to `false` if one wants the algorithm to print details about the iterations and termination of the solver

    - Default is `true`, i.e. by default, there is no output. If one wants to print those information afert solving, the [`print_cnls_model`](@ref) method 
    can be called

* `max_iter::Int` 

    - Maximum number of iterations allowed

    - Default is `100`

* `scaling::Bool` 

    - Set to `true` if one wants the algorithm to work with a constraints jacobian matrix whose rows are scaled (i.e. all constraints gradients vectors are scaled)

    - Default is `false`

* `time_limit::T`

    - Maximum elapsed time (i.e. wall time)

    - Default is `1000`


Tolerances:

* `abs_tol::T`

    - Absolute tolerance for small residuals

    - Default is `eps(T)`

* `rel_tol::T`

    - Relative tolerance used to measure first order criticality and consistency

    - Default is `sqrt(eps(T))`

* `c_tol::T`

    - Tolerance used to measure feasability of the constraints

    - Default is `sqrt(eps(T))`

* `x_tol::T`

    - Tolerance used to measure the distance between two consecutive iterates

    - Default is `sqrt(eps(T))`
"""
function solve!(model::CnlsModel{T}; silent::Bool=true, max_iter::Int=100, scaling::Bool=false, time_limit::T=1e3,
    abs_tol::T=eps(T), rel_tol::T=√abs_tol, c_tol::T=rel_tol, x_tol::T=rel_tol) where {T}

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
    exit_code, x_opt, f_opt, enlsip_info = enlsip(model.starting_point, residuals_evalfunc, constraints_evalfunc, model.nb_parameters, model.nb_residuals, model.nb_eqcons, nb_constraints;
    scaling=scaling, MAX_ITER=max_iter, TIME_LIMIT=time_limit, ε_rel = rel_tol, ε_x = x_tol, ε_c = c_tol, ε_rank=√eps(T))

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