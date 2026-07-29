# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`Enlsip.jl` is a Julia port of the Fortran77 ENLSIP algorithm for solving constrained
nonlinear least-squares problems: minimize `½‖r(x)‖²` subject to nonlinear equality and
inequality constraints plus box bounds on `x`. The method is a Gauss–Newton scheme with an
active-set ("working set") strategy for constraints (Lindström & Wedin, 1988).

## Commands

Run from the package root. The Julia package manager `Pkg` provides the tooling.

- **Run the full test suite:**
  `julia --project -e 'using Pkg; Pkg.test()'`
- **Run tests faster in an existing REPL** (skips re-resolving the sandbox env):
  `julia --project=test test/runtests.jl` after `Pkg.develop`-ing the package into `test/`,
  or simply `include("test/runtests.jl")` from a REPL where `Enlsip` is loaded.
- **Run a single test file** from a REPL:
  `using Enlsip, Test; include("test/problems/HS65.jl")`
  Each test file (`test/internal/*.jl`, `test/problems/*.jl`) is a self-contained `@testset`
  and can be included independently — this is the fastest way to iterate on one area.
- **Build the docs:** `julia --project=docs docs/make.jl`

CI (`.github/workflows/CI.yml`) runs `Pkg.test()` on Julia 1.8 and nightly (Linux). The
minimum supported Julia is 1.8 (`Project.toml` `[compat]`).

## Public API

Users interact with only two things — everything else is internal:

1. **`CnlsModel`** (`src/cnls_model.jl`) — the problem model. Constructed via the keyword
   constructor at `cnls_model.jl:345`. Jacobians are optional; when omitted they are computed
   by `ForwardDiff` (see the one-arg `ResidualsFunction`/`ConstraintsFunction` constructors in
   `src/structures.jl`).
2. **`solve!(model; ...)`** (`src/solver.jl:62`) — runs the solver and mutates `model` in
   place (fills `status_code`, `sol`, `obj_value`, `model_info`). Key kwargs: `silent`,
   `max_iter`, `scaling`, `time_limit`, and the tolerances `abs_tol`/`rel_tol`/`c_tol`/`x_tol`.

Result accessors: `status(model)`, `solution(model)`, `sum_sq_residuals(model)`, plus the
various `*_constraints_values` / `nb_*` query functions exported in `src/cnls_model.jl`.
`status` maps an integer `status_code` through `dict_status_codes` to a `Symbol`
(`:unsolved`, `:found_first_order_stationary_point`, `:failed`,
`:maximum_iterations_exceeded`, `:time_limit_exceeded`).

## Architecture

The module (`src/Enlsip.jl`) includes four files **in order**, and that order is the
dependency chain — later files rely on types/functions from earlier ones:

1. **`structures.jl`** — Core data types. `EvaluationFunction` subtypes
   (`ResidualsFunction`, `ConstraintsFunction`) wrap a value function + Jacobian function and
   count evaluations. `Iteration` / `AbstractIteration` holds all per-iteration state (point,
   search direction, working set, Lagrange multipliers, pseudo-ranks, step length, etc.).
   `WorkingSet` tracks the active-set bookkeeping.

2. **`cnls_model.jl`** — The user-facing `CnlsModel` struct, its constructor, constraint
   builders (`instantiate_constraints_w_bounds` / `_wo_bounds` — bounds are internally folded
   into ordinary inequality constraints), status handling, and result accessors.

3. **`enlsip_functions.jl`** — The solver itself (~2900 lines); this is where the real
   algorithm lives. The entry point is **`enlsip(...)`** (`enlsip_functions.jl:2638`), which
   runs the main iteration loop. The numerical machinery is decomposed into single-purpose
   functions, roughly grouped as:
   - Search direction: `sub_search_direction`, `gn_search_direction` (Gauss–Newton),
     `newton_search_direction`, `search_direction_analys`, `determine_solving_dim`,
     `choose_subspace_dimensions`, `pseudo_rank`.
   - Working set / constraints: `update_working_set`, `init_working_set`,
     `evaluate_violated_constraints`, `check_constraint_deletion`, Lagrange-multiplier
     estimates (`first_/second_lagrange_mult_estimate!`, `minmax_lagrangian_mult`).
   - Line search & step length: `linesearch_constrained`, `compute_steplength`,
     `goldstein_armijo_step`, `upper_bound_steplength`, and the polynomial-root helpers
     (`minrm`, `minrn`, `newton_raphson`, `parameters_rm`) that use `Polynomials.jl`.
   - Penalty weights: `penalty_weight_update`, `euclidean_norm_weight_update`,
     `max_norm_weight_update!`, `min_norm_w!`.
   - Convergence/restart: `check_termination_criteria`, `evaluation_restart!`,
     `check_derivatives`, plus the `print_*` diagnosis helpers.

4. **`solver.jl`** — Thin orchestration layer: `solve!` translates a `CnlsModel` into the raw
   arguments `enlsip` expects, calls it, converts the exit code (`convert_exit_code`), and
   writes results back into the model. Also defines `print_cnls_model`.

The `CnlsModel` → `enlsip` boundary is the key seam: `solve!` is the only place that bridges
the user-friendly model and the low-level solver signature.

## Conventions

- The codebase is parametric on the float type `T<:AbstractFloat` throughout; preserve this
  when editing (e.g. tolerances default to `eps(T)` / `√eps(T)`).
- Unicode identifiers matching the mathematical notation are used freely (`α`, `λ`, `δ`, `ε`,
  `ψ`). Match the surrounding style.
- Mutating functions follow the Julia `!` convention and write results into preallocated
  output arrays passed as arguments (e.g. `res_eval!`, `jacres_eval!`).
- Internal comments are being migrated to English; keep new comments in English.
