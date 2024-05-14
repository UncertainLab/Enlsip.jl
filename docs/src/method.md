# [Method](@id Method)

In this section, a brief description of the method implemented in `Enlsip.jl` is given.

We recall that the problem to solve is of the form:

```math
\begin{aligned}
\min_{x \in \mathbb{R}^n} \quad &  \dfrac{1}{2} \|r(x)\|^2 \\
\text{s.t.} \quad & c_i(x) = 0, \quad i \in \mathcal{E} \\
& c_i(x) \geq 0, \quad i \in \mathcal{I}, \\
\end{aligned}
```

with:

* the residuals $r_i:\mathbb{R}^n\rightarrow\mathbb{R}$ and the constraints $c_i:\mathbb{R}^n\rightarrow\mathbb{R}$ assumed to be $\mathcal{C}^2$ functions;
* norm $\|\cdot\|$ denotes the Euclidean norm.

We introduce the Lagrangian associated to the problem:

```math
L(x,\lambda) = \dfrac{1}{2} \|r(x)\|^2 - \sum_{i \in \mathcal{E} \cup \mathcal{I}} \lambda_i c_i(x),
```

where $\lambda = \left(\lambda_i\right)_{i\in \mathcal{E} \cup \mathcal{I}}$ denotes the vector of Lagrange multipliers.

## Conduct of the algorithm

Starting from a point $x_0$, the algorithm builds a sequence $(x_k)_k$ converging to a first order critical point of the problem.

At a given iteration $k$, a search direction $p_k\in\mathbb{R}^n$ and a steplength $\alpha_k\in[0,1]$ are computed such that the next step is $x_{k+1}:=x_k+\alpha_kp_k$.

An estimate $\lambda_{k}\in\mathbb{R}^{|\mathcal{E}| + |\mathcal{I}|}$ of the vector of Lagrange multipliers is then computed and convergence is tested.

### Search direction

The search direction is the solution of a subproblem derived from an approximation of the original problem.

First, the residuals are linearized in a small neighborhood of current point $x_k$:

$$r(x_k+p)\approx J_kp+r_k,$$

where $J_k$ denotes the Jacobian matrix of the residuals function evaluated at $x_k$ and $r_k:=r(x_k)$ . The resulting objective function of the subproblem corresponds to the Gauss-Newton approximation.

The constraints of the subproblem are formed after a subset of the constraints of the original problem. It contains all the equality constraints but only the inequality constraints that are believed to be satisfied with equality at the solution. This subset, often denoted as the working set in the literature, is updated at every iteration and can be seen as a guess of the optimal active set. For more details on the active set strategy implemented in `Enlsip.jl`, see chapters 5 and 6 of *Practical Optimization*[^GMW82].

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

### Lagrange multipliers estimate

Components of $\lambda_k$ associated to the inequality constraints not in the working set are fixed to zero.

The remaining multipliers, associated to the constraints in the working set, are estimated by the minimum norm solution of the KKT system obtained after cancelling the Lagrange multipliers corresponding to inactive inequality constraints (under strict complementarity). This results into computing the vector $\hat{\lambda}_k^{LS}$ such that:

```math
\hat{\lambda}_k^{LS} \in \arg\min_v \|\hat{A}_k^T v - J_k^Tr_k\|^2.
```

Since this process does not depend from the computation of the search direction, the method implemented in `Enlsip.jl` can be qualified as a primal method.

## [Stopping criteria](@id Stopping)

Termination criteria implemented in `Enlsip.jl` follow the conventions presented in chapter 8 of *Practical Optimization*[^GMW82].

From the Lagrange multipliers estimates $\lambda_k$, we define:

* scalar $\sigma_{min}$: smallest Lagrange multiplier estimate corresponding to an inequality constraint in the working set;
* scalar $\lambda_{max}$: Lagrange multiplier estimate of largest absolute value (among multipliers associated with equality and inequality constraints).

The small positive constants `c_tol`, `rel_tol`, `x_tol` and `abs_tol` are user-specified tolerances. Default values used in `Enlsip.jl` are given in the [Usage](@ref) page (see function [`solve!`](@ref)).

To qualify as a candidate solution, current iterate must first meet the following necessary conditions:

* *Strict Complementarity*: Lagrange multipliers estimates of inequality constraints in the working set are stricly positive;
* *Primal feasability*: $\|\hat{c}_k\| \leq$ `c_tol` and values of the inequality constraints not in the working set are strictly positive at $x_k$;
* *First order criticality*: $\|\nabla_x L(x_k,\lambda_k)\| \leq$ `√rel_tol` $* \left(1+\|J^T_kr_k\|\right)$;
* *Consistency*: $\sigma_{min} \geq$ `rel_tol`$* \lambda_{max}$.

If all criteria above are satisfied, the algorithm stops if one of the following conditions is met:

* *Small distance between the last two setps*: $\|x_{k+1}-x_k\|<$ `x_tol` $ *\|x_k\|$;
* *Small residuals*: $\|r_k\|^2 \leq$ `abs_tol`.

### Convergence

To our knowledge, there is no formal proof of convergence of the method described above, though local convergence with a linear rate is to be expected from the Gauss-Newton paradigm, provided that the initial point is close enough to the solution. Numerical tests confirmed that the efficiency of the method is influenced by the initial point.

[^GMW82]:  P.E. Gill, W. Murray and M.H. Wright, *Practical Optimization*, Academic Press, San Diego, CA, USA, 1982.

[^LW84]: P. Lindström and P.Å. Wedin, *A new linesearch algorithm for nonlinear least squares problems*, Mathematical Programming, vol. 29(3), pages 268-296, 1984.
