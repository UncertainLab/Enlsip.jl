

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