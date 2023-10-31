---
title: 'Enlsip.jl: A Julia optimization package to solve constrained nonlinear least-squares problems'
tags:
  - Julia
  - optimization
  - constrained nonlinear least squares
  - Fortran77
authors:
  - name: Pierre Borie
    orcid: 0009-0000-1043-5057
    affiliation: 1
  - name: Alain Marcotte
    affiliation: 2
  - name: Fabian Bastin
    orcid: 0000-0003-1323-6787
    affiliation: 1
  - name: Stéphane Dellacherie
    affiliation: "2, 3"
affiliations:
 - name: Department of Computer Science and Operations Research, University of Montreal, Montreal, QC, Canada
   index: 1
 - name: Unit of Inflow and Load Forecasting, Hydro-Québec, Montreal, QC, Canada
   index: 2
 - name: Department of Computer Science, UQÀM, Montreal, QC, Canada
   index: 3
date: 31 October 2023
bibliography: paper.bib
---

# Summary

Easy Nonlinear Least Squares Inequality Program (``ENLSIP``[^1]) is the name of an optimization algorithm and an open-source Fortran77 library developed and released in 1988. It implements a nonlinear least squares under nonlinear constraints solver.

This type of problems is mathematically formulated as:
\begin{subequations}\label{eq:cnlls}
         \quad \begin{align}  
                        \quad	\min_{x \in \mathbb{R}^n}        \quad&  \dfrac{1}{2} \sum_{i=1}^{m} r_i(x)^2  \\
                        \text{s.t.}      \quad & c_i(x)=0, \quad 1\leq i \leq q \\
                        & c_i(x) \geq 0, \quad q+1 \leq i \leq \ell,
        \end{align}
\end{subequations}

where the multi-functions $r_i$, often denoted as the residuals, and constraints $c_i$ are two-times differentiable functions. Integers $n,m,q$ and $l$ are the dimensions of the problem.

The ``ENLISP`` solver incorporates a Gauss-Newton type method developed by @lindstromwedin1988. This method also uses an active set strategy to handle the constraints [see @nocedalwright:2006, chapter 16].

[^1]: The source code is available at [https://plato.asu.edu/sub/nonlsq.html](https://plato.asu.edu/sub/nonlsq.html)

# Statement of need

The ``ENLSIP`` Fortran77 library has been successfully used for decades by Hydro-Québec, the main electricity supplier for the province of Quebec in Canada, to calibrate its short-term electricity demand forecast models, which are coded in Fortran90. Since Hydro-Québec is transitioning from Fortran77 to Julia [@Julia:2017] and because its systems are used in a highly critical context, the primary goal of this transition is to ensure that the replacing Julia version reproduces the results given by the original Fortran77 version. The conversion of the above-mentioned ``ENLSIP`` library to Julia is a crucial part of this process.

Nonlinear least squares arise in a variety of model calibration scenarios. Formulation \eqref{eq:cnlls} is particularly relevant in contexts where additional constraints, such as those related to physical models, need to be satisfied. This is due to the high-risk nature of Hydro-Québec's forecasting operations.

Comparison of results and performance on operational Hydro-Québec optimization problems have been conducted using a Julia-Fortran interface and they have shown very good concordance results. We additionaly compared numerical results on nonlinear programming test problems [@hockschittkowski:1980; @lucksanvlcek:1999] to ascertain whether the two versions could yield significantly disparate outcomes or distinct solutions. On the tested problems, we observed no differences in convergence behavior. Furthermore, the obtained solutions did not differ from a predetermined tolerance, the same one we previously employed to ensure the results of our Julia version were consistent with the requirements of Hydro-Québec. This has led us to consider that the current version of our implementation can be published as the Julia package \texttt{Enlsip.jl}.

## From Fortran77 to Julia

Our first motivation to convert ``ENLSIP`` in Julia was to improve reliability, readability and ease of maintenance of the original code. Also, linear algebra tools in Julia, based on [OpenBLAS](http://www.openblas.net), benefit from improved implementations than those of the algorithm by @lindstromwedin1988, based on [MINPACK](https://www.netlib.org/minpack/).
Furthermore, this language is highly convenient for optimization, offering various interface tools such as ``JuMP`` [@JuMP:2017] or ``NLPModels`` [@nlp-models:2020], to model optimization problems. Although these libraries are not currently used in our package, they are under consideration for future developments.
Finally, after conducting some comparison tests on Hydro-Québec operational context problems, we observed that performances of our Julia version were similar to, if not better than, the Fortran77 version in terms of computation time.

## Other nonlinear least-squares packages

Several existing Julia packages can be used to solve nonlinear least-squares problems, such as [NL2sol.jl](https://github.com/macd/NL2sol.jl), [NLS_Solver.jl](https://github.com/vincent-picaud/NLS_Solver.jl) or [CaNNOLeS.jl](https://github.com/JuliaSmoothOptimizers/CaNNOLeS.jl). However, they do not entirely cover the formulation stated in \eqref{eq:cnlls}. Indeed, the first two are designed for unconstrained [@dennisetal:1981] or bound constrained problems and the last one is designed for equality constrained problems [@orbansiquiera:2020].

Although this algorithm may not benefit from state-of-the-art least-squares and nonlinear optimization improvements, its use remains relevant.
Indeed, its application remains very general, covering nonlinearity and non-convexity of the residuals and constraints. Compared to other categories, like the unconstrained case [as discussed in @dennisschnabel:1996, chapter 10], this specific class of least-squares problems with general constraints is, to the best of our knowledge, rarely addressed in the literature.

# Usage

\texttt{Enlsip.jl} can be dowloaded from the Julia package manager by running the following command into the REPL:

```julia
using Pkg 
Pkg.add("Enlsip")
```

Our package provides a basic interface for modeling optimization problems like in \eqref{eq:cnlls}, by passing the residuals, constraints functions and dimensions of the problem.
This is accomplished by creating an instance of our `CnlsModel` structure.
Users can also provide functions to compute Jacobian matrices of residuals and constraints, or they can let the algorithm compute them numerically using forward differences.

We then demonstrate how to use the package by modeling and solving a small dimensions test problem [@hockschittkowski:1980, problem 65].

```julia
using Enlsip

# Dimensions of the problem

n = 3 # number of parameters
m = 3 # number of residuals
l = 1 # number of nonlinear inequality constraints

# Residuals and Jacobian matrix associated
r(x::Vector) = [x[1] - x[2]; (x[1]+x[2]-10.0) / 3.0; x[3]-5.0]

jac_r(x::Vector) = [1. -1. 0; 1/3 1/3 0.; 0. 0. 1.]

# Constraints (one nonlinear inequality and box constraints)
c(x::Vector) = [48.0 - x[1]^2-x[2]^2-x[3]^2]
jac_c(x::Vector) = [ -2x[1] -2x[2] -2x[3]]
x_l = [-4.5, -4.5, -5.0]
x_u = [4.5, 4.5, 5.0] 

# Starting point 
x0 = [-5.0, 5.0, 0.0]

# Instantiate the  model associated with the problem 
model = Enlsip.CnlsModel(r, n, m; jacobian_residuals=jac_r, starting_point=x0,
      ineq_constraints = c, jacobian_ineqcons=jac_c, nb_ineqcons = l,
      x_low=x_l, x_upp=x_u)
```

Once a model has been instantiated, the solver function can be called. If wanted by the user, some details about the iterations can be shown.

```julia
# Call of the `solve!` function
Enlsip.solve!(model)
```

Additional information on how to use the package and examples with test problems from the literature can be found in the [online documentation](https://uncertainlab.github.io/Enlsip.jl/dev)[^2].

[^2]: [https://uncertainlab.github.io/Enlsip.jl](https://uncertainlab.github.io/Enlsip.jl)

# Acknowledgements

The research has been supported by MITACS grants IT25656, IT28724, and IT36208. We would also like to mention that the release of \texttt{Enlsip.jl} results from a close collaboration between the University of Montreal and the Unit of Inflow and Load Forecasting of Hydro-Québec. The work of Fabian Bastin is supported by the Natural Sciences and Engineering Research Council of Canada [Discovery Grant 2022-04400].

# References