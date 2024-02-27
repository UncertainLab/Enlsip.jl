# [Enlsip.jl documentation](@id Home)

## Introduction

Package `Enlsip.jl` is the Julia version of an eponymous Fortran77 library (ENLSIP standing for Easy Nonlinear Least Squares Inequalities Program) designed to solve nonlinear least squares problems under nonlinear constraints.

The optimization method implemented in `Enlsip.jl` was conceived in the late 1980s by two swedish authors named Per Lindström and Per Åke Wedin [^1]. 

It is designed for solve nonlinear least squares problems subject to (s.t.) nonlinear constraints, which can be modeled as the following optimization problem:

```math
\begin{aligned}
\min_{x \in \mathbb{R}^n} \quad &  \dfrac{1}{2} \|r(x)\|^2 \\
\text{s.t.} \quad & c_i(x) = 0, \quad i \in \mathcal{E} \\
& c_i(x) \geq 0, \quad i \in \mathcal{I}, \\
\end{aligned}
```

where:

* the residuals $r_i:\mathbb{R}^n\rightarrow\mathbb{R}$ and the constraints $c_i:\mathbb{R}^n\rightarrow\mathbb{R}$ are assumed to be $\mathcal{C}^2$ functions;
* norm $\|\cdot\|$ denotes the Euclidean norm.

Note that box constraints are modeled as general inequality constraints.

## How to install

To add Enlsip, use Julia's package manager by typing the following command inside the REPL:

```julia
using Pkg
Pkg.add("Enlsip")
```

## How to use

Using `Enlsip.jl` to solve optimization problems consists in, first, instantiating a model and then call the solver on it.

Details and examples with problems from the litterature in the [Usage](@ref) page.

## Description of the algorithm

Starting from a point $x_0$, the algorithm builds a sequence $(x_k)_k$ converging to a first order critical point of the problem.

At a given iteration $k$ , a search direction $p_k\in\mathbb{R}^n$ and a steplength $\alpha_k\in[0,1]$ are computed such that the next step is $x_{k+1}:=x_k+\alpha_kp_k$.

The search direction is the solution of a linear least squares problem under linear constraints.

The objective of the subproblem is obtained after linearizing the residuals in a small neighborhood of current point $x_k$ , resulting into the Gauss-Newton approximation: 
$$r(x_k+p)\approx J_kp+r_k,$$

where $J_k$ denotes the jacobian matrix of the residuals function evaluated at $x_k$ and $r_k:=r(x_k)$ .

The constraints of the subproblem are formed after a subset of the constraints of the original problem, often denoted as the working set in the literature. It contains all the equality constraints but only the inequality constraints that are believed to be satisfied with equality at the solution. This subset, updated at every iteration, can thus be seen as a guess of the optimal active set[^2]. 

Noting $\hat{c}_k$ the vector of constraints in the working set evaluated at $x_k$ and $\hat{A}_k$ its jacobian matrix, the subproblem is then given by:

```math
\begin{aligned}
\min_{p \in \mathbb{R}^n} \quad &  \dfrac{1}{2} \|J_kp+r_k\|^2 \\
\text{s.t.} \quad & \hat{A}_kp + \hat{c}_k= 0.
\end{aligned}
```


[^1]: P. Lindström and P.Å. Wedin, *Gauss-Newton based algorithms for constrained nonlinear least squares problems*, Institute of Information processing, University of Umeå Sweden, 1988.

[^2]: For more details on the active set strategy implemented in `Enlsip.jl`, see chapters 5 and 6 of *Practical Optimization* (P.E. Gill, W. Murray, M.H. Wright, 1982).
