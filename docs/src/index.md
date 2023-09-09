# [Enlsip.jl documentation](@id Home)

Documentation for Enlsip.jl

## Introduction

Package `Enlsip` is the Julia version of an eponymous Fortran77 library (ENLSIP standing for Easy Nonlinear Least Squares Inequalities Program) designed to solve nonlinear least squares problems under nonlinear constraints.

The optimization method implemented in Enlsip was conceived in the late 1980s by two swedish authors named Per Lindström and Per Åke Wedin [^1].

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

## How to use

Using `Enlsip.jl` to solve optimization problems consists in, first, instantiating a model and then call the solver on it. 

Details and examples with problems from the litterature on the [Usage](@ref) and [API](@ref) pages.

## Contents

```@contents
```

[^1]: P. Lindström and P.Å. Wedin, Gauss-Newton based algorithms for constrained nonlinear least squares problems, Institute of Information processing, University of Umeå Sweden, 1988.
