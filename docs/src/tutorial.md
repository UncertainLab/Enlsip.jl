# [Usage](@id Usage)

```@meta
CurrentModule = Enlsip
```

```@setup tutorial
using Enlsip
```

This sections provides details on how to instantiate and solve a constrained least squares problem with `Enlsip.jl`
As a reminder from [Home](@ref), problems to solve are of the following form:

```math
\begin{aligned}
\min_{x \in \mathbb{R}^n} \quad &  \dfrac{1}{2} \|r(x)\|^2 \\
\text{s.t.} \quad & c_i(x) = 0, \quad i \in \mathcal{E} \\
& c_i(x) \geq 0, \quad i \in \mathcal{I}, \\
\end{aligned}
```

with:

* the residuals $r_i:\mathbb{R}^n\rightarrow\mathbb{R}$ and the constraints $c_i:\mathbb{R}^n\rightarrow\mathbb{R}$ assumed to be $\mathcal{C}^2$ functions;
* norm $\|\cdot\|$ denoting the Euclidean norm.
Note that with this formulation, bounds constraints are not distinguished from general inequality constraints. Though, as shown later in this section, they can be provided as vectors of lower and/or upper bounds, which is more convenient for this type of constraints.

It should be borne in mind, however, that the method implemented in `Enlsip.jl` has been conceived for nonlinear problems, as there is no other assumption made about the nature of the residuals and constraints functions, apart from being two-time continously differentiable. The algorithm can still be used to solve linear least squares subject to linear constraints but it will not be as effective as other software where those aspects are taken into account in the design of the optimization method.

## Instantiate a model

Solving a problem with Enlsip is organized in two steps.

First, a model of type [`CnlsModel`](@ref) must be instantiated.

The `CnlsModel` constructor requires the evaluation functions of residuals, constraints, their associated jacobian matrices and dimensions of the problem.

Although the package enables one to create linear unconstrained least squares, it is recommended to use it to solve nonlinear least squares with general constraints.

The three following positional arguments are mandatory to create a model:

* `residuals` : function that computes the vector of residuals
* `nb_parameters` : number of variables
* `nb_residuals` : number of residuals

The following keywords arguments are optional and deal with constraints and Jacobian matrices.
If the Jacobian matrices functions are not provided, they are computed numerically by forward differences using automatic differenciation[^Backend].

[^Backend]: `ForwardDiff.jl` [https://juliadiff.org/ForwardDiff.jl/stable/](https://juliadiff.org/ForwardDiff.jl/stable/)

 Argument             | Details
:---------------------|:----------------------------------------------
 `starting_point`     | initial solution (can be an infeasbile point)
 `jacobian_residuals` | function computing the Jacobian matrix of the residuals
 `eq_constraints`     | function computing the equality constraints
 `jacobian_eqcons`    | function computing the Jacobian matrix of the equality constraints
 `nb_eqcons`          | number of equality constraints
 `ineq_constraints`   | function computing the inequality constraints
 `jacobian_ineqcons`  | function computing the Jacobian matrix of the inequality constraints
 `nb_ineqcons`        | number of inequality constraints
 `x_low`              | vector of lower bounds
 `x_upp`              | vector of upper bounds

It is assumed that the the different functions passed as arguments of the `CnlsModel` constructor are called as `f(x)`, where `x` is a vector of `nb_parameters` elements and `f` is one of the functions `residuals`, `eq_constraints`, `jacobian_eqcons` etc.

## [Solving a model](@id Solving a model)

Then, the `Enlsip` solver can be used by calling the [`solve!`](@ref) function on a instantiated model. By default, the tolerance used in the algorithm is the square root of the relative precision on floating point numbers. For instance, with `Float64`, it will approximately equal `1e-8`.

```@docs
Enlsip.solve!
```

Diagnosis of the conduct of the algorithm can be printed by either setting the `silent` keyword argument of the function [`solve!`](@ref) to `false` or by calling [`print_cnls_model`](@ref) after solving. Here are some details on how to read and understand the different columns of the output:

 Column                        | Description
:------------------------------|:----------------------------------------------
 `iter`                        | iteration number
 `objective`                   | value of the sum of squared residuals (i.e. objective function) at current point
 $\vert\vert$ `active_constraints` $\vert\vert^2$     | value of the sum of squared active constraints at current point
 $\vert\vert$ `p` $\vert\vert$ | norm of the search direction computed at current iteration
 $\alpha$                      | value of the steplength computed at current iteration
 `reduction`                   | reduction in the objective function performed after moving to the next iterate

One can get additional info about termination of the algorithm by calling one of the following functions:

 Name                        |
:----------------------------|
[`solution`](@ref)           |
[`status`](@ref)             |
[`constraints_values`](@ref) |
[`objective_value`](@ref)    |

```@docs
Enlsip.solution
```

```@docs
Enlsip.status
```

```@docs
Enlsip.constraints_values
```

```@docs
Enlsip.objective_value
```

## [Examples](@id Examples)

### Problem 65 from Hock and Schittkowski collection[^HS80]

We show how to implement and solve the following problem:

```math
\begin{aligned}
\min_{x_1, x_2, x_3} \quad & (x_1-x_2)^2 + \dfrac{(x_1+x_2-10)^2}{9}+(x_3-5)^2 \\ 
\text{s.t.} \quad & 48-x_1^2-x_2^2-x_3^2 \geq 0\\
& -4.5\leq x_i \leq 4.5, \quad i=1,2\\
& -5 \leq x_3  \leq 5.
\end{aligned}.
```

The expected optimal solution is $(3.650461821, 3.65046168, 4.6204170507)$.

Associated value of objective function equals $0.9535288567$.

First, we provide the dimensions of the problems.

```@example tutorial
# Dimensions of the problem

n = 3 # number of parameters
m = 3 # number of residuals
nb_eq = 0 # number of equality constraints
nb_constraints = 7 # number of inequality constraints
nothing # hide
```

Then, we define the functions required to compute the residuals, constraints, their respective jacobian matrices and a starting point. In this example, we use the
starting point given in the reference[^HS80], i.e. $(-5, 5, 0)$

```@example tutorial
# Residuals and Jacobian matrix associated
r(x::Vector) = [x[1] - x[2]; (x[1]+x[2]-10.0) / 3.0; x[3]-5.0]

jac_r(x::Vector) = [1. -1. 0;
    1/3 1/3 0.;
    0. 0. 1.]

# Constraints (one equality and box constraints)

c(x::Vector) = [48.0 - x[1]^2-x[2]^2-x[3]^2] # evaluation function for the equality constraint
jac_c(x::Vector) = [ -2x[1] -2x[2] -2x[3]] # Jacobian matrix of the equality constraint

x_l = [-4.5, -4.5, -5.0] # lower bounds
x_u = [4.5, 4.5, 5.0] # upper bounds

# Starting point 
x0 = [-5.0, 5.0, 0.0]
nothing # hide
```

A `CnlsModel` can now be instantiated.

```@example tutorial
# Instantiate a model associated with the problem 
hs65_model = Enlsip.CnlsModel(r, n, m ;starting_point=x0, ineq_constraints = c, 
nb_ineqcons = 1, x_low=x_l, x_upp=x_u, jacobian_residuals=jac_r, jacobian_ineqcons=jac_c)
nothing # hide
```

Finally, the `solve!` function can be called on our model. In this example, keyword arguments remain to default values.

```@example tutorial
Enlsip.solve!(hs65_model)
```

Once `Enlsip` solver has been executed on our problem, a summary of the conduct of the algorithm can be printed by calling [`print_cnls_model`](@ref).

```@example tutorial
Enlsip.print_cnls_model(hs65_model)
```

If one just wants to know about termination of the algorithm, calling [`status`](@ref) will tell if the problem has been successfully solved or not.

```@example tutorial
Enlsip.status(hs65_model)
```

Then, calling [`solution`](@ref) and [`objective_value`](@ref) will respectively return the optimal solution obtained and the value of objective function at that point.

```@example tutorial
hs65_solution = Enlsip.solution(hs65_model)
```

```@example tutorial
hs65_objective = Enlsip.objective_value(hs65_model)
```

The solution obtained is relatively close to the expected optimal solution, although it differs from more than the tolerance used.

```@example tutorial
maximum(abs.(hs65_solution - [3.650461821, 3.65046168, 4.6204170507])) < sqrt(eps(Float64))
```

However, the difference between the objective value obtained with `Enlsip` and the expected one does not exceed the default tolerance.

```@example tutorial
abs(hs65_objective - 0.9535288567) < sqrt(eps(Float64))
```

[^HS80]: W. Hock and K. Schittkowski. *Test Examples for Nonlinear Programming Codes*, volume 187 of Lecture Notes in Economics and Mathematical Systems. Springer, second edition, 1980.
