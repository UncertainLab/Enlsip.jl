# [Usage](@id Usage)

This sections provides details on how to instantiate and solve a constrained least squares problem with `Enlsip.jl`
As a reminder from [Home](@ref), problems to solve are of the following form:

```math
\begin{aligned}
\min \quad &  \dfrac{1}{2} \|r(x)\|^2 \\
\text{s.t.} \quad & c_i(x) = 0, \quad i =1,\ldots,q \\
& c_j(x) \geq 0, \quad j=q+1,\ldots,\ell, \\
\end{aligned}
```

Note that with this formulation, bounds constraints are not distinguished from general inequality constraints. Though, as shown later in this section, they can be provided as vectors of lower and/or upper bounds, which is more convenient for this type of constraints.

Also, the `Enlsip` solver works with double precision float numbers (i.e. type `Float64`).

## Instantiate a model

Solving a problem with Enlsip is organized in two steps.

First, a model must be created by using the [`CnlsModel`](@ref) constructor. This constructor requires the evaluation functions of residuals, constraints, their associated jacobian matrices and dimensions of the problem. 

Although the package enables one to create linear unconstrained least squares, it is recommended to use it to solve nonlinear least squares with general constraints.

The following three arguments are mandatory to create a model:

* `residuals` : function that computes the vector of residuals
    
* `nb_parameters` : number of variables
    
* `nb_residuals` : number of residuals

The following keywords arguments are optionnal and deal with constraints and jacobian matrices computations. If the jacobian matrices functions are not provided, they are computed numerically by forward differences within `Enlsip`.

Argument             | Details
---------------------|----------------------------------------------
`starting_point`     | initial solution (can be an infeasbile point)
`jacobian_residuals` | function computing the jacobian matrix of the residuals
`eq_constraints`     | function computing the equality constraints
`jacobian_eqcons`    | function computing the jacobian matrix of the equality constraints
`nb_eqcons`          | number of equality constraints
`ineq_constraints`   | function computing the inequality constraints
`jacobian_ineqcons`  | function computing the jacobian matrix of the inequality constraints
`nb_ineqcons`        | number of inequality constraints
`x_low`              | vector of lower bounds
`x_upp`              | vector of upper bounds



## Solving a model

Then, the `Enlsip` solver can be used by calling the [`solve!`](@ref) function on a instantiated model. See the [API](@ref) for additionnal information on optionnal arguments that can be passed when calling the solver.


## Examples

### Problem 65 from Hock and Schittkowski collection[^1]

We show how to implement and solve the following problem:

```math
\begin{aligned}
\min \quad & (x_1-x_2)^2 + \dfrac{(x_1+x_2-10)^2}{9}+(x_3-5)^2 \\ 
\text{s.t.} \quad & 48-x_1^2-x_2^2-x_3^2 \geq 0\\
& -4.5\leq x_i \leq 4.5, \quad i=1,2\\
& -5 \leq x_3  \leq 5.
\end{aligned}.
```

First, we provide the dimensions of the problems.

```julia
# Dimensions of the problem

n = 3 # number of parameters
m = 3 # number of residuals
nb_eq = 0 # number of equality constraints
nb_constraints = 7 # number of inequality constraints
```

Then we define the functions required to compute the residuals, constraints, their respective jacobian matrices and a starting point.

```julia
# Residuals and jacobian matrix associated
r(x::Vector) = [x[1] - x[2]; (x[1]+x[2]-10.0) / 3.0; x[3]-5.0]

jac_r(x::Vector) = [1. -1. 0;
    1/3 1/3 0.;
    0. 0. 1.]

# Constraints (one equality and box constraints)

c(x::Vector) = [48.0 - x[1]^2-x[2]^2-x[3]^2] # evaluation function for the equality constraint
jac_c(x::Vector) = [ -2x[1] -2x[2] -2x[3]] # jacobian matrix of the equality constraint

x_l = [-4.5, -4.5, -5.0] # lower bounds
x_u = [4.5, 4.5, 5.0] # upper bounds

# Starting point 
x0 = [-5.0, 5.0, 0.0]
```

Finally, we can instantiate our model and solve it.

```julia
# Instantiate a model associated with the problem 
hs65_model = Enlsip.CnlsModel(r, n, m ;starting_point=x0, ineq_constraints = c, nb_ineqcons = 1, x_low=x_l, x_upp=x_u, 
jacobian_residuals=jac_r, jacobian_ineqcons=jac_c)

# Call of the `solve` function
hs65_sol = Enlsip.solve!(hs65_model,silent=true)
```

[^1]: W. Hock and K. Schittkowski. Test Examples for Nonlinear Programming Codes, volume 187 of Lecture Notes in Economics and Mathematical Systems. Springer, second edition, 1980.
