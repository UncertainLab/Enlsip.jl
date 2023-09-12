# Enlsip.jl

This package is the Julia version of an optimization library originally written in Fortran77 and designed to solve nonlinear least squares problems under nonlinear constraints.
 The optimization method implemented in Enlsip was conceived by two swedish authors, Per Lindstrom and Per Ake Wedin from the Institute of Informatation processing of the University of Umea in Sweden.

Problems that can be solved using Enlsip are modeled as follows:

```math
\begin{aligned}
\min \quad &  \dfrac{1}{2} \|r(x)\|^2 \\
\text{s.t.} \quad & c_i(x) = 0, \quad i =1,\ldots,q \\
& c_j(x) \geq 0, \quad j=q+1,\ldots,\ell, \\
\end{aligned}
```

where $r:\mathbb{R}^n\rightarrow\mathbb{R}^m$, the residuals, and
$c:\mathbb{R}^n\rightarrow\mathbb{R}^{q+\ell}$, concatenation of the constraints, are $\mathcal{C}^1$ multi-functions.

Note that box constraints are modeled as general inequality constraints.

## How to install

To add Enlsip, use Julia's package manager by typing the following command inside the REPL:

```julia
   using Pkg
   Pkg.add("Enlsip")
```

## How to Use Enlsip

Solving a problem with Enlsip is organized in two steps.

First, you need to define to create a model of your problem with the `CnlsModel` structure.

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


### Example with problem 65 from  Hock Schittkowski collection

```julia
# Import Enlsip
using Enlsip

# Dimensions of the problem

n = 3 # number of parameters
m = 3 # number of residuals
nb_eq = 0 # number of equality constraints
nb_constraints = 7 # number of inequality constraints

# Residuals and jacobian matrix associated
r(x::Vector) = [x[1] - x[2]; (x[1]+x[2]-10.0) / 3.0; x[3]-5.0]

jac_r(x::Vector) = [1. -1. 0;
    1/3 1/3 0.;
    0. 0. 1.]

# Constraints (one equality and box constraints)
c(x::Vector) = [48.0 - x[1]^2-x[2]^2-x[3]^2]
jac_c(x::Vector) = [ -2x[1] -2x[2] -2x[3]]
x_l = [-4.5, -4.5, -5.0]
x_u = [4.5, 4.5, 5.0] 

# Starting point 
x0 = [-5.0, 5.0, 0.0]

# Instantiate a model associated with the problem 
hs65_model = Enlsip.CnlsModel(r, n, m ;jacobian_residuals=jac_r, starting_point=x0, ineq_constraints = c, jacobian_ineqcons=jac_c, nb_ineqcons = 1, x_low=x_l, x_upp=x_u)


# Call of the `solve!` function
Enlsip.solve!(hs65_model)

# Print solution and objective value

println("Algorithm successfully terminated : ", Enlsip.status(hs65_model))
println("Optimal solution : ", Enlsip.solution(hs65_model))
println("Optimal objective value : ", ENlsip.objective_value(hs65_model))
```
