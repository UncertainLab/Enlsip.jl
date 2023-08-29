# How to use

Details on how to nstantiate and solve an optmization problem of the following form (as described in [Home](@ref)):

```math
\begin{aligned}
\min \quad &  \dfrac{1}{2} \|r(x)\|^2 \\
\text{s.t.} \quad & c_i(x) = 0, \quad i =1,\ldots,q \\
& c_j(x) \geq 0, \quad j=q+1,\ldots,\ell, \\
\end{aligned}
```

## Instantiate a model

Solving a problem with Enlsip is organized in two steps.

First, a model must be created by using the [`CnlsModel`](@ref) structure.

## Solving a model

Then, once a model is instantiated, the [`solve`](@ref) function may be called.

This function returns an object of type [`CnlsResult`](@ref).

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
hs65_sol = Enlsip.solve(hs65_model,silent=true)
```

[^1]: W. Hock and K. Schittkowski. Test Examples for Nonlinear Programming Codes, volume 187 of Lecture Notes in Economics and Mathematical Systems. Springer, second edition, 1980.
