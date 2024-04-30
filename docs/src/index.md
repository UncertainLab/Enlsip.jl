# [Enlsip.jl documentation](@id Home)

## Introduction

Package `Enlsip.jl` is the Julia version of an eponymous Fortran77 library (ENLSIP standing for Easy Nonlinear Least Squares Inequalities Program) designed to solve nonlinear least squares problems under nonlinear constraints.

The optimization method implemented in `Enlsip.jl` was conceived in the early 1980s by two Swedish authors named Per Lindström and Per Åke Wedin [^LW88].

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

Details and examples with problems from the literature can be found in the [Usage](@ref) page.

## Description of the algorithm

Starting from a point $x_0$, the algorithm builds a sequence $(x_k)_k$ converging to a first order critical point of the problem.

At a given iteration $k$ , a search direction $p_k\in\mathbb{R}^n$ and a steplength $\alpha_k\in[0,1]$ are computed such that the next step is $x_{k+1}:=x_k+\alpha_kp_k$.

### Search direction

The search direction is the solution of a subproblem derived from an approximation of the original problem.

First, the residuals are linearized in a small neighborhood of current point $x_k$:

$$r(x_k+p)\approx J_kp+r_k,$$

where $J_k$ denotes the Jacobian matrix of the residuals function evaluated at $x_k$ and $r_k:=r(x_k)$ . The resulting objective function of the subproblem corresponds to the Gauss-Newton approximation.

The constraints of the subproblem are formed after a subset of the constraints of the original problem. It contains all the equality constraints but only the inequality constraints that are believed to be satisfied with equality at the solution. This subset, often denoted as the working set in the literature, is updated at every iteration and can be seen as a guess of the optimal active set[^2].

Noting $\hat{c}_k$ the vector of constraints in the working set evaluated at $x_k$ and $\hat{A}_k$ its Jacobian matrix, the subproblem is then given by:

```math
\begin{aligned}
\min_{p \in \mathbb{R}^n} \quad &  \dfrac{1}{2} \|J_kp+r_k\|^2 \\
\text{s.t.} \quad & \hat{A}_kp + \hat{c}_k= 0.
\end{aligned}
```

The subproblem is then solved by a null-space type method.

### Steplength

The steplength aims to maintain feasibility  of all of the constraints and to reduce the value of the objective function. In `Enlsip.jl`, this process involves an $\ell_2$-type merit function:

$$\psi_2(x, \mu_k) = \dfrac{1}{2} \|r(x)\|^2 +  \mu_k \sum_{i\in\mathcal{E}} c_i(x)^2+  \mu_k \sum_{i\in\mathcal{I}}  \min(0,c_i(x))^2,$$

where the scalar $\mu_k > 0$ is a penalty parameter updated at every iteration.

Steplength computation is performed by applying a linesearch method on function $\psi_2$. This consists in minimizing the merit function along the direction from $x_k$ to $x_k+p_k$ , i.e. finding $\alpha_k\in[0,1]$ such that

$$\alpha_k \in \arg\min_{\alpha \in [0,1]} \psi_2(x_k+\alpha_kp_k, \mu_k).$$

The authors of the Fortran77 version of ENLSIP developed a linesearch method in which an approximate, but acceptable, minimizer of the merit function is computed[^LW84].

### Convergence

To our knowledge, there is no formal proof of convergence of the method described above, though local convergence with a linear rate should be expected from the Gauss-Newton paradigm, provided that the initial point is close enough to the solution and that the optimal active does not change. Numerical tests confirmed that the efficiency of the method is influenced by the initial point.

## Bug reports and contributions

As this package is a conversion from Fortran77 to Julia, there might be some bugs that we did not encountered yet, so if you think you found one, you can open an [issue](https://github.com/UncertainLab/Enlsip.jl/issues) to report it.

Issues can also be opened to discuss about eventual suggestions of improvement.

[^LW88]: P. Lindström and P.Å. Wedin, *Gauss-Newton based algorithms for constrained nonlinear least squares problems*, Institute of Information processing, University of Umeå Sweden, 1988.

[^2]: For more details on the active set strategy implemented in `Enlsip.jl`, see chapters 5 and 6 of *Practical Optimization* (P.E. Gill, W. Murray, M.H. Wright, 1982).

[^LW84]: P. Lindström and P.Å. Wedin, *A new linesearch algorithm for nonlinear least squares problems*, Mathematical Programming, vol. 29(3), pages 268-296, 1984.
