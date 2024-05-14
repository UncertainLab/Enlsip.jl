# [Enlsip.jl documentation](@id Home)

## Introduction

Package `Enlsip.jl` is the Julia version of an eponymous Fortran77 library (ENLSIP standing for Easy Nonlinear Least Squares Inequalities Program) designed to solve nonlinear least-squares problems under nonlinear constraints.

The optimization method implemented in `Enlsip.jl` was conceived in the early 1980s by two Swedish authors named Per Lindström and Per Åke Wedin [^LW88].

It is designed for solve nonlinear least-squares problems subject to (s.t.) nonlinear constraints, which can be modeled as the following optimization problem:

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

`Enlsip.jl` incorporates an iterative method computing a first order critical point of the problem. A brief description of the method and the stopping criteria is given in [Method](@ref).

## Bug reports and contributions

As this package is a conversion from Fortran77 to Julia, there might be some bugs that we did not encountered yet, so if you think you found one, you can open an [issue](https://github.com/UncertainLab/Enlsip.jl/issues) to report it.

Issues can also be opened to discuss about eventual suggestions of improvement.

[^LW88]: P. Lindström and P.Å. Wedin, *Gauss-Newton based algorithms for constrained nonlinear least squares problems*  , Institute of Information processing, University of Umeå Sweden, 1988.
