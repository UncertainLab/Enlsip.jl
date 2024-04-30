# Enlsip.jl

[![](https://img.shields.io/badge/docs-stable-green.svg)](https://uncertainlab.github.io/Enlsip.jl/stable/) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://uncertainlab.github.io/Enlsip.jl/dev/) 

Package `Enlsip.jl` is the Julia version of ENLSIP, an open source algorithm originally written in Fortran77 and designed to solve nonlinear least-squares problems subject to nonlinear constraints.
The optimization method implemented in ENLSIP is described in
 
> Per Lindström and Per-Åke Wedin, *Gauss Newton based algorithms for constrained nonlinear least squares problems*.
> Technical Report S-901 87, Institute of Information processing, University of Umeå, Sweden, 1988.

The source code of the Fortran77 library is available at [https://plato.asu.edu/sub/nonlsq.html](https://plato.asu.edu/sub/nonlsq.html).

Problems that can be solved using `Enlsip.jl` are modeled as follows:

```math
\begin{aligned}
\min_{x \in \mathbb{R}^n} \quad &  \dfrac{1}{2} \|r(x)\|^2 \\
\text{s.t.} \quad & c_i(x) = 0, \quad i \in \mathcal{E} \\
& c_i(x) \geq 0, \quad i \in \mathcal{I}, \\
& l_i \leq x_i \leq u_i, \quad i \in \{1,\ldots,n\}, 
\end{aligned}
```

where:

* the residuals $r_i:\mathbb{R}^n\rightarrow\mathbb{R}$ and the constraints $c_i:\mathbb{R}^n\rightarrow\mathbb{R}$ are assumed to be $\mathcal{C}^2$ functions;
* norm $\|\|\cdot\|\|$ denotes the Euclidean norm;
* $l$ and $u$ are respectively vectors of lower and upper bounds.

In the formulation above, bound constraints are written seperately but they are treated as classical inequality constraints in the method implemented in ENLSIP.

## How to install

To add Enlsip, use Julia's package manager by typing the following command inside the REPL:

```julia
   using Pkg
   Pkg.add("Enlsip")
```

## How to Use

Solving a problem with Enlsip is organized in two steps.

First, you need to create a model of your problem with the `CnlsModel` structure.

### Creating a model

An object of type `CnlsModel` can be created using a constructor, whose arguments are the following:

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

### Solving a model

Then, once your model is instantiated, you can call the `solve!` function to solve your problem.

### Example with problem 65 from Hock Schittkowski collection[^HS80]

```julia
# Import Enlsip
using Enlsip

# Dimensions of the problem

n = 3 # number of parameters
m = 3 # number of residuals

# Residuals and jacobian matrix associated
r(x::Vector) = [x[1] - x[2]; (x[1]+x[2]-10.0) / 3.0; x[3]-5.0]

jac_r(x::Vector) = [1. -1. 0;
    1/3 1/3 0.;
    0. 0. 1.]

# Constraints (one nonlinear inequality and box constraints)
c(x::Vector) = [48.0 - x[1]^2-x[2]^2-x[3]^2]
jac_c(x::Vector) = [ -2x[1] -2x[2] -2x[3]]
x_l = [-4.5, -4.5, -5.0]
x_u = [4.5, 4.5, 5.0] 

# Starting point 
x0 = [-5.0, 5.0, 0.0]

# Instantiate a model associated with the problem 
hs65_model = Enlsip.CnlsModel(r, n, m ;jacobian_residuals=jac_r, starting_point=x0,
ineq_constraints = c, jacobian_ineqcons=jac_c, nb_ineqcons = 1, x_low=x_l, x_upp=x_u)


# Call of the `solve!` function
Enlsip.solve!(hs65_model)

# Print solution and objective value

println("Algorithm termination status: ", Enlsip.status(hs65_model))
println("Optimal solution: ", Enlsip.solution(hs65_model))
println("Optimal objective value: ", Enlsip.objective_value(hs65_model))
```
[^HS80]: W. Hock and K. Schittkowski. *Test Examples for Nonlinear Programming Codes*, volume 187 of Lecture Notes in Economics and Mathematical Systems. Springer, second edition, 1980.
