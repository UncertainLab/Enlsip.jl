export AbstractCnlsModel, CnlsModel
export total_nb_constraints, status, solution, objective_value, dict_status_codes

#= Structures where the first two fields are the functions evaluating residuals or constraints and the associated jacobian matrix
=#


abstract type EvaluationFunction end

mutable struct ResidualsFunction <: EvaluationFunction 
    reseval
    jacres_eval
    nb_reseval::Int64
    nb_jacres_eval::Int64
end


function ResidualsFunction(eval, jac_eval)
    ResidualsFunction(eval, jac_eval, 0, 0)
end

function ResidualsFunction(eval)
    num_jac_eval(x::Vector) = jac_forward_diff(eval,x)
    ResidualsFunction(eval, num_jac_eval,0,0)
end

mutable struct ConstraintsFunction <: EvaluationFunction 
    conseval
    jaccons_eval
    nb_conseval::Int64
    nb_jaccons_eval::Int64
end

function ConstraintsFunction(eval, jac_eval)
    ConstraintsFunction(eval, jac_eval, 0, 0)
end

function ConstraintsFunction(eval)
    num_jac_eval(x::Vector) = jac_forward_diff(eval,x)
    ConstraintsFunction(eval, num_jac_eval,0,0)
end

#= Functions to compute in place the residuals, constraints and jacobian matrices of a given EvaluationFunction =#

function res_eval!(r::ResidualsFunction, x::Vector, rx::Vector)
    rx[:] = r.reseval(x)
    r.nb_reseval += 1
    return
end

function jacres_eval!(r::ResidualsFunction, x::Vector, J::Matrix)
    J[:] = r.jacres_eval(x)
    r.nb_jacres_eval += 1
    return
end

function cons_eval!(c::ConstraintsFunction, x::Vector, cx::Vector)
    cx[:] = c.conseval(x)
    c.nb_conseval += 1
    return
end

function jaccons_eval!(c::ConstraintsFunction, x::Vector, A::Matrix)
    A[:] = c.jaccons_eval(x)
    c.nb_jaccons_eval += 1
    return
end

# Auxialiry funcion to define EvaluationFunction with numerical jacobian
function jac_forward_diff(h_eval, x::Vector)

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

"""
    AbstractCnlsModel

Abstract type for [`CnlsModel`](@ref) structure.
"""
abstract type AbstractCnlsModel end


struct CNLSModel <: AbstractCnlsModel
    residuals::ResidualsFunction
    constraints::ConstraintsFunction
    starting_point::Vector
    nb_parameters::Int64
    nb_residuals::Int64
    nb_eqcons::Int64
    nb_cons::Int64
end



function CNLSModel(
    residuals=nothing,
    nb_parameters::Int64=0,
    nb_residuals::Int64=0;
    starting_point::Vector=zeros(Float64, nb_parameters),
    jacobian_residuals=nothing,
    eq_constraints=nothing,
    jacobian_eqcons=nothing,
    nb_eqcons::Int64=0,
    ineq_constraints=nothing,
    jacobian_ineqcons=nothing,
    nb_ineqcons::Int64=0,
    x_low::Vector=fill!(Vector{Float64}(undef,nb_parameters), -Inf),
    x_upp::Vector=fill!(Vector{Float64}(undef,nb_parameters), Inf))
    

    # Assertion test on residuals
    @assert(typeof(residuals) <: Function, "The argument res_func must be a function")
    @assert(nb_parameters > 0 && nb_residuals > 0, "The number of parameters and number of residuals must be strictly positive")
    
    # Assertion tests on constraints
    @assert(eq_constraints !== nothing || ineq_constraints !== nothing || any(isfinite,x_low) || any(isfinite,x_upp), "There must be at least one constraint")
    @assert(!(eq_constraints===nothing && nb_eqcons != 0) || !(typeof(eq_constraints <: Function) && nb_eqcons == 0), "Incoherent definition of equality constraints")
    @assert(!(ineq_constraints===nothing && nb_ineqcons != 0) || !(typeof(ineq_constraints <: Function) && nb_ineqcons == 0), "Incoherent definition of inequality constraints")

    residuals_evalfunc = (jacobian_residuals === nothing ? ResidualsFunction(residuals) : ResidualsFunction(residuals, jacobian_residuals))

    if all(!isfinite,vcat(x_low,x_upp))
        constraints_evalfunc = instantiate_constraints_wo_bounds(eq_constraints, jacobian_eqcons, ineq_constraints, jacobian_ineqcons)
    else
        constraints_evalfunc = instantiate_constraints_w_bounds(eq_constraints, jacobian_eqcons, ineq_constraints, jacobian_ineqcons, x_low, x_upp)
    end

    nb_constraints = nb_eqcons + nb_ineqcons + count(isfinite, x_low) + count(isfinite,x_upp)

    return CNLSModel(residuals_evalfunc, constraints_evalfunc, starting_point, nb_parameters, nb_residuals, nb_eqcons, nb_constraints)
end

#= TODO

Faire une structure mutable CnlsModel qui contient toutes les infos du modèle avec possibilité de les ajouter par l'utilisateur
    Doit y avoir les fonctions pour ajouter ces infos comme sur Jump (sauf que c'est limité à des fonctions qui retournent des vecteurs)

Faire une autre structure non mutable EnlsipModel, qui remplace le CnlsModel AsbtractCnlsResult

Faire de CnlsResult des attributs du nouveau CnlsModel
=#


"""
    CnlsModel

Structure modeling an instance of a constrainted nonlinear least squares problem.



    * `residuals` : function that computes the vector of residuals
    
    * `nb_parameters` : number of variables
    
    * `nb_residuals` : number of residuals

    * `stating_point` : initial solution
    
    * `jacobian_residuals` : function that computes the jacobian matrix of the residuals
    
    * `eq_constraints` : function that computes the vector of equality constraints
    
    * `jacobian_eqcons` : function that computes the jacobian matrix of the equality constraints
    
    * `nb_eqcons` : number of equality constraints
    
    * `ineq_constraints` : function that computes the vector of inequality constraints
    
    * `jacobian_ineqcons` : function that computes the jacobian matrix of the inequality constraints
    
    * `nb_ineqcons` : number of inequality constraints
    
    * `x_low` and `x_upp` : respectively vectors of lower and upper bounds

    * `status_code` : integer indicating the solving status of the model. Call function [`status`](@ref)
"""
mutable struct CnlsModel <: AbstractCnlsModel
    residuals
    nb_parameters::Int
    nb_residuals::Int
    starting_point::Vector
    jacobian_residuals
    eq_constraints
    jacobian_eqcons
    nb_eqcons::Int
    ineq_constraints
    jacobian_ineqcons
    nb_ineqcons::Int
    x_low::Vector
    x_upp::Vector
    status_code::Int
    sol::Vector
    obj_value::Float64
end

function convert_exit_code(code::Int64)
    
    status_code = 0
    if code > 0
        status_code = 1
    elseif code == -2
        status_code = code
    else
        status_code = -1
    end

    return status_code
end

const dict_status_codes = Dict(
    0 => :unsolved,
    1 => :successfully_solved,
    -1 => :failed,
    -2 => :maximum_iterations_exceeded,
)


"""
    status(model)

Returns readable information on the solving status of `model`. 

See also:[`CnlsModel`](@ref).
"""
status(model::CnlsModel) = dict_status_codes[model.status_code]

"""
    solution(model)

Once the given `model` has been solved, this function returns the optimal solution, or last solution obtained if no convergence, as a `Vector` of approriate dimension.

See also: [`solve!`](@ref) and [`CnlsModel`](@ref).
"""
solution(model::CnlsModel) = model.sol

"""
    objective_value(model)

Once the given `model` has been solved, returns the value of the objective function, i.e. sum of squared residuals functions, computed at the optimal solution.
If no convergence, this value is computed at the last solution obtained.

See also: [`CnlsModel`](@ref) and [`solution`](@ref).
"""
objective_value(model::CnlsModel) = model.obj_value

"""
    total_nb_constraints(model)

Returns the total number of constraints, i.e. equalities, inequalities and bounds, of the given `model`.

See also: [`CnlsModel`](@ref).
"""
total_nb_constraints(model::CnlsModel) = model.nb_eqcons + model.nb_ineqcons + count(isfinite, model.x_low) + count(isfinite, model.x_upp)

"""
    model = CnlsModel(residuals, nb_parameters, nb_residuals)

Constructor for [`CnlsModel`](@ref)

The following arguments are mandatory to instantiate a model:

    * `residuals` : function that computes the vector of residuals
    
    * `nb_parameters` : number of variables
    
    * `nb_residuals` : number of residuals
    

The following can be provided as optionnal arguments:

    * `stating_point` : initial solution (default is a vector of zeros of appropriate dimension)
    
    * `jacobian_residuals` : function that computes the jacobian matrix of the residuals. If not passed as argument, it is computed numericcaly by forward differences
    
    * `eq_constraints` : function that computes the vector of equality constraints
    
    * `jacobian_eqcons` : function that computes the jacobian matrix of the equality constraints. If not passed as argument, it is computed numericcaly by forward differences
    
    * `nb_eqcons` : number of equality constraints
    
    * `ineq_constraints` : function that computes the vector of inequality constraints
    
    * `jacobian_ineqcons` : function that computes the jacobian matrix of the inequality constraints. If not passed as argument, it is computed numericcaly by forward differences
    
    * `nb_ineqcons` : number of inequality constraints
    
    * `x_low` and `x_upp` : respectively vectors of lower and upper bounds
"""
function CnlsModel(
    residuals=nothing,
    nb_parameters::Int64=0,
    nb_residuals::Int64=0;
    starting_point::Vector=zeros(Float64, nb_parameters),
    jacobian_residuals=nothing,
    eq_constraints=nothing,
    jacobian_eqcons=nothing,
    nb_eqcons::Int64=0,
    ineq_constraints=nothing,
    jacobian_ineqcons=nothing,
    nb_ineqcons::Int64=0,
    x_low::Vector=fill!(Vector{Float64}(undef,nb_parameters), -Inf),
    x_upp::Vector=fill!(Vector{Float64}(undef,nb_parameters), Inf))
    

    # Assertion test on residuals
    @assert(typeof(residuals) <: Function, "The residuals argument must be a function")
    @assert(nb_parameters > 0 && nb_residuals > 0, "The number of parameters and number of residuals must be strictly positive")
    
    # Assertion tests on constraints
    @assert(eq_constraints !== nothing || ineq_constraints !== nothing || any(isfinite,x_low) || any(isfinite,x_upp), "There must be at least one constraint")
    @assert(!(eq_constraints===nothing && nb_eqcons != 0) || !(typeof(eq_constraints <: Function) && nb_eqcons == 0), "Incoherent definition of equality constraints")
    @assert(!(ineq_constraints===nothing && nb_ineqcons != 0) || !(typeof(ineq_constraints <: Function) && nb_ineqcons == 0), "Incoherent definition of inequality constraints")

    rx0 = residuals(starting_point)
    initial_obj_value = dot(rx0,rx0)
    return CnlsModel(residuals, nb_parameters, nb_residuals, starting_point, jacobian_residuals, 
    eq_constraints, jacobian_eqcons, nb_eqcons, ineq_constraints, jacobian_ineqcons, nb_ineqcons, x_low, x_upp,
    0, starting_point, initial_obj_value)
end




function box_constraints(x_low::Vector, x_upp::Vector)

    n = size(x_low,1)
    @assert(n == size(x_upp,1),"Bounds vectors must have same length")

    no_x_low = all(!isfinite,x_low)
    no_x_upp = all(!isfinite,x_upp)

    @assert(!(no_x_low && no_x_upp), "Bounds vectors are assumed to contain at least one finite element")

    if no_x_low && !no_x_upp
        cons_w_ubounds(x::Vector) = filter(isfinite,x_upp-x)
        jaccons_w_ubounds(x::Vector) = Matrix{eltype(x_upp)}(-I,n,n)[filter(i-> isfinite(x_upp[i]),1:n),:]
        return cons_w_ubounds, jaccons_w_ubounds
    
    elseif !no_x_low && no_x_upp
        cons_w_lbounds(x::Vector) = filter(isfinite, x-x_low)
        jaccons_w_lbounds(x::Vector) = Matrix{eltype(x_low)}(I,n,n)[filter(i-> isfinite(x_low[i]),1:n),:]
        return cons_w_lbounds, jaccons_w_lbounds
    
    else
        cons_w_bounds(x::Vector) = vcat(filter(isfinite, x-x_low), filter(isfinite, x_upp-x))
        jaccons_w_bounds(x::Vector) = vcat(Matrix{eltype(x_low)}(I,n,n)[filter(i-> isfinite(x_low[i]),1:n),:], Matrix{eltype(x_upp)}(-I,n,n)[filter(i-> isfinite(x_upp[i]),1:n),:])
        return cons_w_bounds, jaccons_w_bounds
    end
end

function instantiate_constraints_w_bounds(eq_constraints, jacobian_eqcons, ineq_constraints, jacobian_ineqcons, x_low, x_upp)

    bounds_func, jac_bounds_func = box_constraints(x_low, x_upp)

    if eq_constraints !== nothing && ineq_constraints !== nothing
        cons(x::Vector) = vcat(eq_constraints(x), ineq_constraints(x), bounds_func(x))

        if jacobian_eqcons !== nothing && jacobian_ineqcons !== nothing
            jac_cons(x::Vector) = vcat(jacobian_eqcons(x), jacobian_ineqcons(x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_cons)

        elseif jacobian_eqcons !== nothing && jacobian_ineqcons === nothing
            jac_cons_ineqnum(x::Vector) = vcat(jacobian_eqcons(x), jac_forward_diff(ineq_constraints, x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_ineqnum)

        elseif jacobian_eqcons === nothing && jacobian_ineqcons !== nothing
            jac_cons_eqnum(x) = vcat(jac_forward_diff(eq_constraints,x), jacobian_ineqcons(x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_eqnum)
        else
            jac_consnum(x::Vector) =  vcat(jac_forward_diff(eq_constraints,x), jac_forward_diff(ineq_constraints, x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_consnum)
        end

    elseif eq_constraints !== nothing && ineq_constraints === nothing
        eq_cons(x::Vector) = vcat(eq_constraints(x), bounds_func(x))
        
        if jacobian_eqcons === nothing
            jac_eqconsnum(x::Vector) = vcat(jac_forward_diff(eq_constraints,x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(eq_cons, jac_eqconsnum)
        else
            jac_eqcons(x::Vector) = vcat(jacobian_eq_cons(x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(eq_cons, jac_eqcons)
        end

    elseif eq_constraints === nothing && ineq_constraints !== nothing
        ineq_cons(x::Vector) = vcat(ineq_constraints(x), bounds_func(x))
        
        if jacobian_ineqcons === nothing
            jac_ineqconsnum(x::Vector) = vcat(jac_forward_diff(ineq_constraints,x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(ineq_cons, jac_ineqconsnum)
        else
            jac_ineqcons(x::Vector) = vcat(jacobian_ineqcons(x), jac_bounds_func(x))
            constraints_evalfunc = ConstraintsFunction(ineq_cons, jac_ineqcons)
        end
    
    else
        constraints_evalfunc = ConstraintsFunction(bounds_func, jac_bounds_func)
    end
    return constraints_evalfunc
end

function instantiate_constraints_wo_bounds(eq_constraints, jacobian_eqcons, ineq_constraints, jacobian_ineqcons)

    if eq_constraints !== nothing && ineq_constraints !== nothing
        cons(x::Vector) = vcat(eq_constraints(x), ineq_constraints(x))

        if jacobian_eqcons !== nothing && jacobian_ineqcons !== nothing
            jac_cons(x::Vector) = vcat(jacobian_eqcons(x), jacobian_ineqcons(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_cons)

        elseif jacobian_eqcons !== nothing && jacobian_ineqcons === nothing
            jac_cons_ineqnum(x::Vector) = vcat(jacobian_eqcons(x), jac_forward_diff(ineq_constraints, x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_ineqnum)

        elseif jacobian_eqcons === nothing && jacobian_ineqcons !== nothing
            jac_cons_eqnum(x) = vcat(jac_forward_diff(eq_constraints,x), jacobian_ineqcons(x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_eqnum)
        else
            jac_consnum(x::Vector) =  vcat(jac_forward_diff(eq_constraints,x), jac_forward_diff(ineq_constraints, x))
            constraints_evalfunc = ConstraintsFunction(cons, jac_consnum)
        end

    elseif eq_constraints !== nothing && ineq_constraints === nothing
        constraints_evalfunc = (jacobian_eqcons === nothing ? ConstraintsFunction(eq_constraints) : ConstraintsFunction(eq_constraints, jacobian_eqcons))
        
    elseif eq_constraints === nothing && ineq_constraints !== nothing
        constraints_evalfunc = (jacobian_ineqcons === nothing ? ConstraintsFunction(ineq_constraints) : ConstraintsFunction(ineq_constraints, jacobian_ineqcons))
    end

    return constraints_evalfunc
end