export AbstractEnlsipModel, EnlsipModel

#=

    Iteration

Summarizes the useful informations about an iteration of the algorithm

* `x` : Departure point of the iteration 

* `p` : Descent direction

* `rx` : vector of size `m`, contains value of residuals at `x` 

* `cx` : vector of size `l`, contains value of constraints at `x`

* `t` : Number of constraints in current working set (ie constraints considered active)

* `α` : Value of steplength

* `λ` : Vector of size `t`, containts Lagrange multipliers estimates

* `rankA` : pseudo rank of matrix `A`, jacobian of active constraints

* `rankJ2` : pseudo rank of matrix `J2`, block extracted from `J`, jacobian of residuals

* `b_gn` : right handside of the linear system solved to compute first part of `p`

* `d_gn` :  right handside of the linear system solved to compute second part of `p`

* `predicted_reduction` : predicted linear progress

* `progress` :  reduction in the objective function

* `β` : scalar used to estimate convergence factor

* `restart` : indicate if current iteration is a restart step or no

* `first` : indicate if current iteration is the first one or no

* `add` : indicate if a constraint has been added to the working set 

* `del` : indicate if a constraint has been deleted from the working set

* `index_del` : index of the constraint that has been deleted from working set (`0` if no deletion)

* `code` : Its value caracterizes the method used to compute the search direction `p`

    - `1` represents Gauss-Newton method

    - `-1` represents Subspace minimization

    - `2`  represents Newton method

* `nb_newton_steps` : number of search direction computed using the method of Newton
=#
mutable struct Iteration
    x::Vector
    p::Vector
    rx::Vector
    cx::Vector
    t::Int64
    α::Float64
    index_α_upp::Int64
    λ::Vector
    w::Vector
    rankA::Int64
    rankJ2::Int64
    dimA::Int64
    dimJ2::Int64
    b_gn::Vector
    d_gn::Vector
    predicted_reduction::Float64
    progress::Float64
    grad_res::Float64
    speed::Float64
    β::Float64
    restart::Bool
    first::Bool
    add::Bool
    del::Bool
    index_del::Int64
    code::Int64
    nb_newton_steps::Int64
end


Base.copy(s::Iteration) = Iteration(s.x, s.p, s.rx, s.cx, s.t, s.α, s.index_α_upp, s.λ, s.w, s.rankA, s.rankJ2, s.dimA, s.dimJ2, s.b_gn, s.d_gn, 
s.predicted_reduction, s.progress, s.grad_res, s.speed, s.β, s.restart, s.first, s.add, s.del, s.index_del, s.code, s.nb_newton_steps)

#=
    Constraint

Struct used to represent the active constraints

Fields are the useful informations about active constraints at a point x :

* `cx` : Vector of size t, contains values of constraints in current working set

* `A` : Matrix of size `t` x `t`, jacobian matrix of constraints in current working set

* `scaling` : Boolean indicating if internal scaling of `cx` and `A` is done 

* `diag_scale` : Vector of size `t`, contains the diagonal elements of the scaling matrix if internal scaling is done 

    - The i-th element equals ``\\dfrac{1}{\\|\\nabla c_i(x)\\|}`` for ``i = 1,...,t``, which is the inverse of the length of `A` i-th row 
    - Otherwise, it contains the length of each row in the matrix `A`
=#
mutable struct Constraint
    cx::Vector{Float64}
    A::Matrix{Float64}
    scaling::Bool
    diag_scale::Vector{Float64}
end

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

# EVSCAL 
# Scale jacobian matrix of active constraints A and active constraints evaluation vector cx if so indicated (ie if scale different from 0) by forming vectors :
# diag*A and diag*cx
# where diag is an array of dimension whose i-th element equals either ||∇c_i(x)|| or  (1/||∇c_i(x)|) depending on wether scaling is done or not. 
# The vectors are formed by modifying in place matrix A and vector cx 

function evaluate_scaling!(C::Constraint)

    t = size(C.A, 1)
    ε_rel = eps(eltype(C.cx))
    C.diag_scale = zeros(t)
    for i = 1:t
        row_i = norm(C.A[i, :])
        C.diag_scale[i] = row_i
        if C.scaling
            if abs(row_i) < ε_rel
                row_i = 1.0
            end
            C.A[i, :] /= row_i
            C.cx[i] /= row_i
            C.diag_scale[i] = 1.0 / row_i
        end
    end
    return
end


#=
    WorkingSet

In ENLSIP, the working-set is a prediction of the set of active constraints at the solution

It is updated at every iteration thanks to a Lagrangian multipliers estimation

Fields of this structure summarize infos about the qualification of the constraints, i.e. :

* `q` : number of equality constraints

* `t` : number of constraints in current working set (all equalities and some inequalities considered to be active at the solution)

* `l` : total number of constraints (i.e. equalities and inequalities)

* active :

    - `Vector` of size `l`

    - first `t` elements are indeces of constraints in working set sorted in increasing order, other elements equal `0`

* inactive : 

    - `Vector` of size `l-q`

    - first `l-t` elements are indeces of constraints not in working set sorted in increasing order, other elements equal `0`

=#
mutable struct WorkingSet
    q::Int64
    t::Int64
    l::Int64
    active::Vector{Int64}
    inactive::Vector{Int64}
end

# Equivalent Fortran : DELETE in dblreduns.f
# Moves the active constraint number s to the inactive set

function delete_constraint!(W::WorkingSet, s::Int64)

    l, t = W.l, W.t

    # Ajout de la contrainte à l'ensemble inactif
    W.inactive[l-t+1] = W.active[s]
    sort!(@view W.inactive[1:l-t+1])

    # Réorganisation de l'ensemble actif
    for i = s:t-1
        W.active[i] = W.active[i+1]
    end
    W.active[t] = 0
    W.t -= 1
    return
end

# Equivalent Fortran : ADDIT in dblreduns.f
# Add the inactive constraint nulber s to the active s

function add_constraint!(W::WorkingSet, s::Int64)

    l, t = W.l, W.t
    # s-th inactive constraint moved from inactive to active set
    W.active[t+1] = W.inactive[s]
    sort!(@view W.active[1:t+1])
    # Inactive set reorganized
    for i = s:l-t-1
        W.inactive[i] = W.inactive[i+1]
    end
    W.inactive[l-t] = 0
    W.t += 1
    return
end

abstract type AbstractEnlsipModel end


"""
    EnlsipModel

Structure modeling an instance of a constrainted nonlinear least squares problem 
"""
struct EnlsipModel <: AbstractEnlsipModel
    residuals::ResidualsFunction
    constraints::ConstraintsFunction
    starting_point::Vector
    nb_parameters::Int64
    nb_residuals::Int64
    nb_eqcons::Int64
    nb_cons::Int64
end


"""
    model = EnlsipModel(residuals,nb_parameters,nb_residuals;
    starting_point,jacobian_residuals=nothing,eq_constraints=nothing,jacobian_eqcons=nothing,nb_eqcons,ineq_constraints=nothing,jacobian_ineqcons=nothing,nb_ineqcons,x_low,x_upp)

Constructor for [`EnlsipModel`](@ref)

Arguments are the following:

    * `residuals` : function that computes the vector of residuals
    
    * `nb_parameters` : number of variables
    
    * `nb_residuals` : number of residuals
    
    * `stating_point` : initial solution
    
    * `jacobian_residuals` : function that computes the jacobian matrix of the residuals. If not passed as argument, it is computed numericcaly by forward differences
    
    * `eq_constraints` : function that computes the vector of equality constraints
    
    * `jacobian_eqcons` : function that computes the jacobian matrix of the equality constraints. If not passed as argument, it is computed numericcaly by forward differences
    
    * `nb_eqcons` : number of equality constraints
    
    * `ineq_constraints` : function that computes the vector of inequality constraints
    
    * `jacobian_ineqcons` : function that computes the jacobian matrix of the inequality constraints. If not passed as argument, it is computed numericcaly by forward differences
    
    * `nb_ineqcons` : number of inequality constraints
    
    * `x_low` and `x_upp` : respectively vectors of lower and upper bounds
"""
function EnlsipModel(
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

    return EnlsipModel(residuals_evalfunc, constraints_evalfunc, starting_point, nb_parameters, nb_residuals, nb_eqcons, nb_constraints)
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


    
        
    
