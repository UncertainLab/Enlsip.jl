export AbstractCnlsModel, CnlsModel
export constraints_values, equality_constraints_values, inequality_constraints_values, bounds_constraints_values, total_nb_constraints, nb_equality_constraints, nb_inequality_constraints, nb_lower_bounds, nb_upper_bounds
export status, solution, sum_sq_residuals, dict_status_codes

#= Structures where the first two fields are the functions evaluating residuals or constraints and the associated jacobian matrix
=#


abstract type EvaluationFunction end

mutable struct ResidualsFunction <: EvaluationFunction 
    reseval
    jacres_eval
    nb_reseval::Int
    nb_jacres_eval::Int
end


function ResidualsFunction(eval, jac_eval)
    ResidualsFunction(eval, jac_eval, 0, 0)
end

# function ResidualsFunction(eval)
#     num_jac_eval(x::Vector{<:AbstractFloat}) = jac_forward_diff(eval,x)
#     ResidualsFunction(eval, num_jac_eval,0,0)
# end

ResidualsFunction(eval) = ResidualsFunction(eval, x::Vector -> ForwardDiff.jacobian(eval, x))

mutable struct ConstraintsFunction <: EvaluationFunction 
    conseval
    jaccons_eval
    nb_conseval::Int
    nb_jaccons_eval::Int
end

function ConstraintsFunction(eval, jac_eval)
    ConstraintsFunction(eval, jac_eval, 0, 0)
end

# function ConstraintsFunction(eval)
#     num_jac_eval(x::Vector{<:AbstractFloat}) = jac_forward_diff(eval,x)
#     ConstraintsFunction(eval, num_jac_eval,0,0)
# end

ConstraintsFunction(eval) = ConstraintsFunction(eval, x::Vector -> ForwardDiff.jacobian(eval, x))

#= Functions to compute in place the residuals, constraints and jacobian matrices of a given EvaluationFunction =#

function res_eval!(r::ResidualsFunction, x::Vector{T}, rx::Vector{T}) where {T<:AbstractFloat}
    rx[:] = r.reseval(x)
    r.nb_reseval += 1
    return
end

function jacres_eval!(r::ResidualsFunction, x::Vector{T}, J::Matrix{T}) where {T<:AbstractFloat}
    J[:] = r.jacres_eval(x)
    r.nb_jacres_eval += 1
    return
end

function cons_eval!(c::ConstraintsFunction, x::Vector{T}, cx::Vector{T}) where {T<:AbstractFloat}
    cx[:] = c.conseval(x)
    c.nb_conseval += 1
    return
end

function jaccons_eval!(c::ConstraintsFunction, x::Vector{T}, A::Matrix{T}) where {T<:AbstractFloat}
    A[:] = c.jaccons_eval(x)
    c.nb_jaccons_eval += 1
    return
end

# Auxialiry funcion to define EvaluationFunction with numerical jacobian
function jac_forward_diff(h_eval, x::Vector{<:AbstractFloat})

    δ = sqrt(eps(eltype(x)))
    hx = h_eval(x)
    n = size(x,1)
    m = size(hx,1)

    Jh = zeros(eltype(hx),(m,n))

    for j=1:n
        δ_j = max(abs(x[j]), 1.0) * δ
        e_j = [(i == j ? 1.0 : 0.0) for i = 1:n]
        x_forward = x + δ_j * e_j
        hx_forward = h_eval(x_forward)
        Jh[:,j] = (hx_forward - hx) / δ_j
    end
    return Jh
end

#=
    ExecutionInfo

Summarizes information that are given after termination of the algorithm

* iterations_detail: List of `DisplayedInfo`, each element of the list contains the details about a given iteration. The list is stored in chronological order (1st element -> 1st iterration, last element -> last_iteration)

* nb_function_evaluations : number of residuals and constraints evaluations performed during the execution of the algorithm

* nb_jacobian_evaluations : same as above but for Jacobian matrices evaluations

* solving_time : time of execution in seconds
=#
struct ExecutionInfo
    iterations_detail::Vector{DisplayedInfo}
    nb_function_evaluations::Int
    nb_jacobian_evaluations::Int
    solving_time::Float64
end

ExecutionInfo() = ExecutionInfo([DisplayedInfo()], 0, 0, 0.0)
"""
    AbstractCnlsModel

Abstract type for [`CnlsModel`](@ref) structure.
"""
abstract type AbstractCnlsModel{T<:AbstractFloat} end

"""
    CnlsModel{T} where {T<:AbstractFloat}

Structure modeling an instance of a constrainted nonlinear least squares problem.

This structure contains the following attributes:

    * `residuals` : function that computes the vector of residuals
    
    * `nb_parameters::Int` : number of variables
    
    * `nb_residuals::Int` : number of residuals

    * `stating_point::Vector{T}` : initial solution
    
    * `jacobian_residuals` : function that computes the jacobian matrix of the residuals
    
    * `eq_constraints` : function that computes the vector of equality constraints
    
    * `jacobian_eqcons` : function that computes the jacobian matrix of the equality constraints
    
    * `nb_eqcons` : number of equality constraints
    
    * `ineq_constraints` : function that computes the vector of inequality constraints
    
    * `jacobian_ineqcons` : function that computes the jacobian matrix of the inequality constraints
    
    * `nb_ineqcons::Int` : number of inequality constraints
    
    * `x_low::Vector{T}` and `x_upp::Vector{T}` : respectively vectors of lower and upper bounds

    * `status_code::Int` : integer indicating the solving status of the model
"""
mutable struct CnlsModel{T} <: AbstractCnlsModel{T}
    residuals
    nb_parameters::Int
    nb_residuals::Int
    starting_point::Vector{T}
    jacobian_residuals
    eq_constraints
    jacobian_eqcons
    nb_eqcons::Int
    ineq_constraints
    jacobian_ineqcons
    nb_ineqcons::Int
    x_low::Vector{T}
    x_upp::Vector{T}
    constraints_scaling::Bool
    status_code::Int
    sol::Vector{T}
    obj_value::T
    model_info::ExecutionInfo
end

function convert_exit_code(code::Int)
    
    status_code = 0
    if code > 0
        status_code = 1
    elseif code == -2 || code == -11
        status_code = code
    else
        status_code = -1
    end

    return status_code
end

const dict_status_codes = Dict(
    0 => :unsolved,
    1 => :found_first_order_stationary_point,
    -1 => :failed,
    -2 => :maximum_iterations_exceeded,
    -11 => :time_limit_exceeded
)


"""
    status(model)

This functions returns a `Symbol` that gives brief information on the solving status of `model`.

If a model has been instantiated but the solver has not been called yet, it will return `:unsolved`.

Once the solver has been called and if a first order stationary point satisfying the convergence criteria has been computed, it will return `:found_first_order_stationary_point`.

If the algorithm met an abnormall termination criteria, it will return one of the following:

* `:failed` : the algorithm encoutered a numerical error that triggered termination

* `:maximum_iterations_exceeded` : a solution could not be reached within the maximum number of iterations

* `:time_limit_exceeded` : the algorithm stopped because solving time exceeded the time limit
"""
status(model::CnlsModel) = dict_status_codes[model.status_code]

"""
    solution(model)

Once the given `model` has been solved, this function returns the optimal solution, or last solution obtained if no convergence, as a `Vector` of approriate dimension.
"""
solution(model::CnlsModel) = model.sol

"""
    sum_sq_residuals(model)

Once the given `model` has been solved, returns the value of the objective function, i.e. sum of squared residuals functions, computed at the optimal solution.
If no convergence, this value is computed at the last solution obtained.
"""
sum_sq_residuals(model::CnlsModel) = model.obj_value

nb_equality_constraints(model::CnlsModel) = model.nb_eqcons

nb_inequality_constraints(model::CnlsModel) = model.nb_ineqcons

nb_lower_bounds(model::CnlsModel) = count(isfinite, model.x_low)

nb_upper_bounds(model::CnlsModel) = count(isfinite, model.x_upp)

"""
    total_nb_constraints(model)

Returns the total number of constraints, i.e. equalities, inequalities and bounds, of the given `model`.

See also: [`CnlsModel`](@ref).
"""
total_nb_constraints(model::CnlsModel) = nb_equality_constraints(model) + nb_inequality_constraints(model) + nb_lower_bounds(model) + nb_upper_bounds(model)


"""
    equality_constraints_values(model)

Returns the vector of equality constraints values at the solution of `model` (if they are any).
"""
function equality_constraints_values(model::CnlsModel)
    sol = solution(model)
    hx = Vector{eltype(sol)}(undef, nb_equality_constraints(model))
    if model.eq_constraints !== nothing
        hx[:] = model.eq_constraints(sol)
    end
    return hx
end

"""
    inequality_constraints_values(model)

Returns the vector of inequality constraints values at the solution of `model` (if they are any).
"""
function inequality_constraints_values(model::CnlsModel)
    sol = solution(model)
    gx = Vector{eltype(sol)}(undef, nb_inequality_constraints(model))
    if model.ineq_constraints !== nothing
        gx[:] = model.ineq_constraints(sol)
    end
    return gx
end

"""
    bounds_constraints_values(model)

Returns the vector of box constraints values at the solution `xₛ` of `model` (if they are any).   

If `xₗ` and `xᵤ` are respectively the vectors of lower and upper bounds, it will return `[xₛ-xₗ; xᵤ-xₛ]`.
"""
function bounds_constraints_values(model::CnlsModel)
    sol = solution(model)
    return vcat(sol-model.x_low, model.x_upp-sol)
end

"""
    constraints_values(model)

Computes values of all the constraints in `model` at the solution. 

The vector returned is the concatenation of equalities, inequalities and box constraints (in that order).

For instance, let `xₛ` be the solution found. If functions `h`, `g` compute equality and inequality constraints and `xₗ, xᵤ` are vectors of lower and lower bounds,
it will return `[h(xₛ); g(xₛ); xₛ-xₗ; xᵤ-xₛ]`.

If one wants to compute each type of constraints seperately, see [`equality_constraints_values`](@ref), [`inequality_constraints_values`](@ref) and [`bounds_constraints_values`](@ref).
"""
function constraints_values(model::CnlsModel)
    sol = solution(model)
    q, ℓ = model.nb_eqcons, model.nb_ineqcons

    cx = Vector{eltype(sol)}(undef, q+ℓ+2*length(sol))

    if model.nb_eqcons > 0
        cx[1:q] = equality_constraints_values(model)
    end
    if model.nb_ineqcons > 0
        cx[q+1:q+ℓ] = inequality_constraints_values(model)
    end
    if nb_lower_bounds(model) + nb_upper_bounds(model) > 0
        cx[q+ℓ+1:end] = bounds_constraints_values(model)
    end
    return cx
end

"""
    model = CnlsModel(residuals, nb_parameters, nb_residuals)

Constructor for [`CnlsModel`](@ref).

* Positional arguments

    - `residuals` : function that computes the vector of residuals
    
    - `nb_parameters` : number of variables
    
    - `nb_residuals` : number of residuals
    

* Keywords arguments :

    - `stating_point::Vector{T}` : initial solution (default is a vector of zeros of appropriate dimension)
    
    - `jacobian_residuals` : function that computes the jacobian matrix of the residuals. If not passed as argument, it is computed numericcaly by forward differences
    
    - `eq_constraints` : function that computes the vector of equality constraints
    
    - `jacobian_eqcons` : function that computes the jacobian matrix of the equality constraints. If not passed as argument, it is computed numericcaly by forward differences
    
    - `nb_eqcons::Int` : number of equality constraints
    
    - `ineq_constraints` : function that computes the vector of inequality constraints
    
    - `jacobian_ineqcons` : function that computes the jacobian matrix of the inequality constraints. If not passed as argument, it is computed numericcaly by forward differences
    
    - `nb_ineqcons::Int` : number of inequality constraints
    
    - `x_low::Vector{T}` and `x_upp::Vector{T}` : respectively vectors of lower and upper bounds
"""
function CnlsModel(
    residuals=nothing,
    nb_parameters::Int=0,
    nb_residuals::Int=0;
    starting_point::Vector{T}=zeros(nb_parameters),
    jacobian_residuals=nothing,
    eq_constraints=nothing,
    jacobian_eqcons=nothing,
    nb_eqcons::Int=0,
    ineq_constraints=nothing,
    jacobian_ineqcons=nothing,
    scaling::Bool=false,
    nb_ineqcons::Int=0,
    x_low::Vector{T}=fill!(Vector{eltype(starting_point)}(undef,nb_parameters), -Inf),
    x_upp::Vector{T}=fill!(Vector{eltype(starting_point)}(undef,nb_parameters), Inf)) where {T}
    

    # Assertion test on residuals
    @assert(typeof(residuals) <: Function, "A function evaluating residuals must be provided")
    @assert(nb_parameters > 0 && nb_residuals > 0, "The number of parameters and number of residuals must be strictly positive")
    
    # Assertion tests on constraints
    @assert(eq_constraints !== nothing || ineq_constraints !== nothing || any(isfinite,x_low) || any(isfinite,x_upp), "There must be at least one constraint")
    @assert(!(eq_constraints===nothing && nb_eqcons != 0) || !(typeof(eq_constraints <: Function) && nb_eqcons == 0), "Incoherent definition of equality constraints")
    @assert(!(ineq_constraints===nothing && nb_ineqcons != 0) || !(typeof(ineq_constraints <: Function) && nb_ineqcons == 0), "Incoherent definition of inequality constraints")

    rx0 = residuals(starting_point)
    initial_obj_value = dot(rx0,rx0)

    initial_info = ExecutionInfo()
    return CnlsModel(residuals, nb_parameters, nb_residuals, starting_point, jacobian_residuals, 
    eq_constraints, jacobian_eqcons, nb_eqcons, ineq_constraints, jacobian_ineqcons, nb_ineqcons, x_low, x_upp,
    scaling, 0, starting_point, initial_obj_value, initial_info)
end


function box_constraints(x_low::Vector{T}, x_upp::Vector{T}) where {T<:AbstractFloat}

    n = size(x_low,1)
    @assert(n == size(x_upp,1),"Bounds vectors must have same length")

    no_x_low = all(!isfinite,x_low)
    no_x_upp = all(!isfinite,x_upp)

    @assert(!(no_x_low && no_x_upp), "Bounds vectors are assumed to contain at least one finite element")

    if no_x_low && !no_x_upp
        cons_w_ubounds(x::Vector{<:AbstractFloat}) = filter(isfinite,x_upp-x)
        jaccons_w_ubounds(x::Vector{<:AbstractFloat}) = Matrix{eltype(x_upp)}(-I,n,n)[filter(i-> isfinite(x_upp[i]),1:n),:]
        return cons_w_ubounds, jaccons_w_ubounds
    
    elseif !no_x_low && no_x_upp
        cons_w_lbounds(x::Vector{<:AbstractFloat}) = filter(isfinite, x-x_low)
        jaccons_w_lbounds(x::Vector{<:AbstractFloat}) = Matrix{eltype(x_low)}(I,n,n)[filter(i-> isfinite(x_low[i]),1:n),:]
        return cons_w_lbounds, jaccons_w_lbounds
    
    else
        cons_w_bounds(x::Vector{<:AbstractFloat}) = vcat(filter(isfinite, x-x_low), filter(isfinite, x_upp-x))
        jaccons_w_bounds(x::Vector{<:AbstractFloat}) = vcat(Matrix{eltype(x_low)}(I,n,n)[filter(i-> isfinite(x_low[i]),1:n),:], Matrix{eltype(x_upp)}(-I,n,n)[filter(i-> isfinite(x_upp[i]),1:n),:])
        return cons_w_bounds, jaccons_w_bounds
    end
end

# Returns a ConstraintsEvaluation function for problems with bound constraints

function instantiate_constraints_w_bounds(eq_constraints, jacobian_eqcons, ineq_constraints, jacobian_ineqcons, x_low, x_upp)

    bounds_func, jac_bounds_func = box_constraints(x_low, x_upp)

    # Equality and inequality constraints in the model
    if eq_constraints !== nothing && ineq_constraints !== nothing
        cons(x::Vector) = vcat(eq_constraints(x), ineq_constraints(x), bounds_func(x))

        # Jacobian provided 
        if jacobian_eqcons !== nothing && jacobian_ineqcons !== nothing
            jac_cons(x::Vector) = vcat(jacobian_eqcons(x), jacobian_ineqcons(x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_cons)
        # Jacobian not provided for inequality constraints
        elseif jacobian_eqcons !== nothing && jacobian_ineqcons === nothing
            ad_jac_ineqcons(x::Vector) = vcat(jacobian_eqcons(x), ForwardDiff.jacobian(ineq_constraints, x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(cons, ad_jac_ineqcons)
        # Jacobian not provided for equality constraints
        elseif jacobian_eqcons === nothing && jacobian_ineqcons !== nothing
            ad_jac_eqcons(x) = vcat(ForwardDiff.jacobian(eq_constraints,x), jacobian_ineqcons(x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_eqnum)
        # Jacobian not provided
        else
            ad_jac_cons(x::Vector) =  vcat(ForwardDiff.jacobian(eq_constraints,x), ForwardDiff.jacobian(ineq_constraints, x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(cons, ad_jac_cons)
        end
    # No inequality constraints
    elseif eq_constraints !== nothing && ineq_constraints === nothing
        eq_cons(x::Vector) = vcat(eq_constraints(x), bounds_func(x))
        # Jacobian not provided
        if jacobian_eqcons === nothing
            jac_eqconsnum(x::Vector) = vcat(ForwardDiff.jacobian(eq_constraints,x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(eq_cons, jac_eqconsnum)
        # Jacobian not provided
        else
            jac_eqcons(x::Vector) = vcat(jacobian_eq_cons(x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(eq_cons, jac_eqcons)
        end
    # No equality constraints
    elseif eq_constraints === nothing && ineq_constraints !== nothing
        ineq_cons(x::Vector) = vcat(ineq_constraints(x), bounds_func(x))
        # Jacobian not provided
        if jacobian_ineqcons === nothing
            jac_ineqconsnum(x::Vector) = vcat(ForwardDiff.jacobian(ineq_constraints,x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(ineq_cons, jac_ineqconsnum)
        # Jacobian provided
        else
            jac_ineqcons(x::Vector) = vcat(jacobian_ineqcons(x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(ineq_cons, jac_ineqcons)
        end
    
    else
        constraints_evalfunc = ConstraintsFunction(bounds_func, jac_bounds_func)
    end
    return constraints_evalfunc
end

# Returns a ConstraintsEvaluation function for problems without bound constraints
function instantiate_constraints_wo_bounds(eq_constraints, jacobian_eqcons, ineq_constraints, jacobian_ineqcons)

    if eq_constraints !== nothing && ineq_constraints !== nothing
        cons(x::Vector) = vcat(eq_constraints(x), ineq_constraints(x))

        if jacobian_eqcons !== nothing && jacobian_ineqcons !== nothing
            jac_cons(x::Vector) = vcat(jacobian_eqcons(x), jacobian_ineqcons(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_cons)

        elseif jacobian_eqcons !== nothing && jacobian_ineqcons === nothing
            jac_cons_ineqnum(x::Vector) = vcat(jacobian_eqcons(x), ForwardDiff.jacobian(ineq_constraints, x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_ineqnum)

        elseif jacobian_eqcons === nothing && jacobian_ineqcons !== nothing
            ad_jac_eqcons(x) = vcat(ForwardDiff.jacobian(eq_constraints,x), jacobian_ineqcons(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_eqnum)
        else
            ad_jac_cons(x::Vector) =  vcat(ForwardDiff.jacobian(eq_constraints,x), ForwardDiff.jacobian(ineq_constraints, x))
            constraints_evalfunc = ConstraintsFunction(cons, ad_jac_cons)
        end

    elseif eq_constraints !== nothing && ineq_constraints === nothing
        constraints_evalfunc = (jacobian_eqcons === nothing ? ConstraintsFunction(eq_constraints) : ConstraintsFunction(eq_constraints, jacobian_eqcons))
        
    elseif eq_constraints === nothing && ineq_constraints !== nothing
        constraints_evalfunc = (jacobian_ineqcons === nothing ? ConstraintsFunction(ineq_constraints) : ConstraintsFunction(ineq_constraints, jacobian_ineqcons))
    end

    return constraints_evalfunc
end
