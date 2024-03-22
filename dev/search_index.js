var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [Enlsip]","category":"page"},{"location":"reference/#Enlsip.AbstractCnlsModel","page":"Reference","title":"Enlsip.AbstractCnlsModel","text":"AbstractCnlsModel\n\nAbstract type for CnlsModel structure.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Enlsip.CnlsModel","page":"Reference","title":"Enlsip.CnlsModel","text":"CnlsModel\n\nStructure modeling an instance of a constrainted nonlinear least squares problem.\n\nThis structure contains the following attributes:\n\n* `residuals` : function that computes the vector of residuals\n\n* `nb_parameters` : number of variables\n\n* `nb_residuals` : number of residuals\n\n* `stating_point` : initial solution\n\n* `jacobian_residuals` : function that computes the jacobian matrix of the residuals\n\n* `eq_constraints` : function that computes the vector of equality constraints\n\n* `jacobian_eqcons` : function that computes the jacobian matrix of the equality constraints\n\n* `nb_eqcons` : number of equality constraints\n\n* `ineq_constraints` : function that computes the vector of inequality constraints\n\n* `jacobian_ineqcons` : function that computes the jacobian matrix of the inequality constraints\n\n* `nb_ineqcons` : number of inequality constraints\n\n* `x_low` and `x_upp` : respectively vectors of lower and upper bounds\n\n* `status_code` : integer indicating the solving status of the model.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Enlsip.CnlsModel-Union{Tuple{}, Tuple{Any}, Tuple{T}, Tuple{Any, Int64}, Tuple{Any, Int64, Int64}} where T","page":"Reference","title":"Enlsip.CnlsModel","text":"model = CnlsModel(residuals, nb_parameters, nb_residuals)\n\nConstructor for CnlsModel.\n\nPositional arguments\nresiduals : function that computes the vector of residuals\nnb_parameters : number of variables\nnb_residuals : number of residuals\nKeywords arguments :\nstating_point : initial solution (default is a vector of zeros of appropriate dimension)\njacobian_residuals : function that computes the jacobian matrix of the residuals. If not passed as argument, it is computed numericcaly by forward differences\neq_constraints : function that computes the vector of equality constraints\njacobian_eqcons : function that computes the jacobian matrix of the equality constraints. If not passed as argument, it is computed numericcaly by forward differences\nnb_eqcons : number of equality constraints\nineq_constraints : function that computes the vector of inequality constraints\njacobian_ineqcons : function that computes the jacobian matrix of the inequality constraints. If not passed as argument, it is computed numericcaly by forward differences\nnb_ineqcons : number of inequality constraints\nx_low and x_upp : respectively vectors of lower and upper bounds\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.bounds_constraints_values-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.bounds_constraints_values","text":"bounds_constraints_values(model)\n\nReturns the vector of box constraints values at the solution xₛ of model (if they are any).   \n\nIf xₗ and xᵤ are respectively the vectors of lower and upper bounds, it will return [xₛ-xₗ; xᵤ-xₛ].\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.constraints_values-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.constraints_values","text":"constraints_values(model)\n\nComputes values of all the constraints in model at the solution. \n\nThe vector returned is the concatenation of equalities, inequalities and box constraints (in that order).\n\nFor instance, let xₛ be the solution found. If functions h, g compute equality and inequality constraints and xₗ, xᵤ are vectors of lower and lower bounds, it will return [h(xₛ); g(xₛ); xₛ-xₗ; xᵤ-xₛ].\n\nIf one wants to compute each type of constraints seperately, see equality_constraints_values, inequality_constraints_values and bounds_constraints_values.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.equality_constraints_values-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.equality_constraints_values","text":"equality_constraints_values(model)\n\nReturns the vector of equality constraints values at the solution of model (if they are any).\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.inequality_constraints_values-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.inequality_constraints_values","text":"inequality_constraints_values(model)\n\nReturns the vector of inequality constraints values at the solution of model (if they are any).\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.objective_value-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.objective_value","text":"objective_value(model)\n\nOnce the given model has been solved, returns the value of the objective function, i.e. sum of squared residuals functions, computed at the optimal solution. If no convergence, this value is computed at the last solution obtained.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.print_cnls_model","page":"Reference","title":"Enlsip.print_cnls_model","text":"print_cnls_model(model,io)\n\nOne can call this function to print information about an instance model (see CnlsModel). \n\nIf model has just been instantiated but not solved, it will print general information about the model, such as the dimensions of the residuals, parameters and constraints.\n\nAfter calling the solve! method, the output will be enhanced with details about the iterations performed during the execution of the algorithm.\n\nThe following info are also printed:\n\nnumber of iterations\ntotal number of function and Jacobian matrix evaluations for both residuals and contraints\nsolving time in seconds\nvalue of the objective function found by the algorithm\ntermination status (see status)\n\n\n\n\n\n","category":"function"},{"location":"reference/#Enlsip.solution-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.solution","text":"solution(model)\n\nOnce the given model has been solved, this function returns the optimal solution, or last solution obtained if no convergence, as a Vector of approriate dimension.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.solve!-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.solve!","text":"solve!(model)\n\nOnce a CnlsModel has been instantiated, this function solves the optimzation problem associated by using the method implemented in the Enlsip solver.\n\nKeywords arguments:\n\nsilent::Bool \nSet to false if one wants the algorithm to print details about the iterations and termination of the solver\nDefault value is true, i.e. by default, there is no output. If one wants to print those information afert solving, the print_cnls_model method \ncan be called.\nmax_iter::Int \nMaximum number of iterations allowed\nDefault value is set to 100\nscaling::Bool \nSet to true if one wants the algorithm to work with a constraints jacobian matrix whose rows are scaled (i.e. all constraints gradients vectors are scaled)\nDefault value is set to false\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.status-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.status","text":"status(model)\n\nThis functions returns a Symbol that gives brief information on the solving status of model.\n\nIf a model has been instantiated but the solver has not been called yet, it will return :unsolved.\n\nOnce the solver has been called and if a first order critical point satisfying the convergence criteria has been computed, it will return :successfully_solved.\n\nIf the algorithm met an abnormall termination criteria, it will return one of the following:\n\n:failed : the algorithm encoutered a numerical error that triggered termination\n:maximum_iterations_exceeded : a solution could not be reached within the maximum number of iterations.\n\n\n\n\n\n","category":"method"},{"location":"reference/#Enlsip.total_nb_constraints-Tuple{CnlsModel}","page":"Reference","title":"Enlsip.total_nb_constraints","text":"total_nb_constraints(model)\n\nReturns the total number of constraints, i.e. equalities, inequalities and bounds, of the given model.\n\nSee also: CnlsModel.\n\n\n\n\n\n","category":"method"},{"location":"#Home","page":"Home","title":"Enlsip.jl documentation","text":"","category":"section"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Package Enlsip.jl is the Julia version of an eponymous Fortran77 library (ENLSIP standing for Easy Nonlinear Least Squares Inequalities Program) designed to solve nonlinear least squares problems under nonlinear constraints.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The optimization method implemented in Enlsip.jl was conceived in the early 1980s by two Swedish authors named Per Lindström and Per Åke Wedin [LW88].","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is designed for solve nonlinear least squares problems subject to (s.t.) nonlinear constraints, which can be modeled as the following optimization problem:","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\nmin_x in mathbbR^n quad   dfrac12 r(x)^2 \ntextst quad  c_i(x) = 0 quad i in mathcalE \n c_i(x) geq 0 quad i in mathcalI \nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"where:","category":"page"},{"location":"","page":"Home","title":"Home","text":"the residuals r_imathbbR^nrightarrowmathbbR and the constraints c_imathbbR^nrightarrowmathbbR are assumed to be mathcalC^2 functions;\nnorm cdot denotes the Euclidean norm.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note that box constraints are modeled as general inequality constraints.","category":"page"},{"location":"#How-to-install","page":"Home","title":"How to install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To add Enlsip, use Julia's package manager by typing the following command inside the REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"Enlsip\")","category":"page"},{"location":"#How-to-use","page":"Home","title":"How to use","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Using Enlsip.jl to solve optimization problems consists in, first, instantiating a model and then call the solver on it.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Details and examples with problems from the literature can be found in the Usage page.","category":"page"},{"location":"#Description-of-the-algorithm","page":"Home","title":"Description of the algorithm","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Starting from a point x_0, the algorithm builds a sequence (x_k)_k converging to a first order critical point of the problem.","category":"page"},{"location":"","page":"Home","title":"Home","text":"At a given iteration k , a search direction p_kinmathbbR^n and a steplength alpha_kin01 are computed such that the next step is x_k+1=x_k+alpha_kp_k.","category":"page"},{"location":"#Search-direction","page":"Home","title":"Search direction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The search direction is the solution of a subproblem derived from an approximation of the original problem.","category":"page"},{"location":"","page":"Home","title":"Home","text":"First, the residuals are linearized in a small neighborhood of current point x_k:","category":"page"},{"location":"","page":"Home","title":"Home","text":"r(x_k+p)approx J_kp+r_k","category":"page"},{"location":"","page":"Home","title":"Home","text":"where J_k denotes the Jacobian matrix of the residuals function evaluated at x_k and r_k=r(x_k) . The resulting objective function of the subproblem corresponds to the Gauss-Newton approximation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The constraints of the subproblem are formed after a subset of the constraints of the original problem. It contains all the equality constraints but only the inequality constraints that are believed to be satisfied with equality at the solution. This subset, often denoted as the working set in the literature, is updated at every iteration and can be seen as a guess of the optimal active set[2].","category":"page"},{"location":"","page":"Home","title":"Home","text":"Noting hatc_k the vector of constraints in the working set evaluated at x_k and hatA_k its Jacobian matrix, the subproblem is then given by:","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\nmin_p in mathbbR^n quad   dfrac12 J_kp+r_k^2 \ntextst quad  hatA_kp + hatc_k= 0\nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"The subproblem is then solved by a null-space type method.","category":"page"},{"location":"#Steplength","page":"Home","title":"Steplength","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The steplength aims to maintain feasibility  of all of the constraints and to reduce the value of the objective function. In Enlsip.jl, this process involves an ell_2-type merit function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"psi_2(x mu_k) = dfrac12 r(x)^2 +  mu_k sum_iinmathcalE c_i(x)^2+  mu_k sum_iinmathcalI  min(0c_i(x))^2","category":"page"},{"location":"","page":"Home","title":"Home","text":"where the scalar mu_k  0 is a penalty parameter updated at every iteration.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Steplength computation is performed by applying a linesearch method on function psi_2. This consists in minimizing the merit function along the direction from x_k to x_k+p_k , i.e. finding alpha_kin01 such that","category":"page"},{"location":"","page":"Home","title":"Home","text":"alpha_k in argmin_alpha in 01 psi_2(x_k+alpha_kp_k mu_k)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The authors of the Fortran77 version of ENLSIP developed a linesearch method in which an approximate, but acceptable, minimizer of the merit function is computed[LW84].","category":"page"},{"location":"#Convergence","page":"Home","title":"Convergence","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To our knowledge, there is no proof of convergence of the method described above, though local convergence with a linear rate should be expected from the Gauss-Newton paradigm, provided that the initial point is close enough to the solution and that the optimal active does not change. Numerical tests confirmed that the efficiency of the method is influenced by the initial point.","category":"page"},{"location":"#Bug-reports-and-contributions","page":"Home","title":"Bug reports and contributions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As this package is a conversion from Fortran77 to Julia, there might be some bugs that we did not encountered yet, so if you think you found one, you can open an issue to report it.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Issues can also be opened to discuss about eventual suggestions of improvement.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[LW88]: P. Lindström and P.Å. Wedin, Gauss-Newton based algorithms for constrained nonlinear least squares problems, Institute of Information processing, University of Umeå Sweden, 1988.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[2]: For more details on the active set strategy implemented in Enlsip.jl, see chapters 5 and 6 of Practical Optimization (P.E. Gill, W. Murray, M.H. Wright, 1982).","category":"page"},{"location":"","page":"Home","title":"Home","text":"[LW84]: P. Lindström and P.Å. Wedin, A new linesearch algorithm for nonlinear least squares problems, Mathematical Programming, vol. 29(3), pages 268-296, 1984.","category":"page"},{"location":"tutorial/#Usage","page":"Usage","title":"Usage","text":"","category":"section"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"CurrentModule = Enlsip","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"using Enlsip","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"This sections provides details on how to instantiate and solve a constrained least squares problem with Enlsip.jl As a reminder from Home, problems to solve are of the following form:","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"beginaligned\nmin_x in mathbbR^n quad   dfrac12 r(x)^2 \ntextst quad  c_i(x) = 0 quad i in mathcalE \n c_i(x) geq 0 quad i in mathcalI \nendaligned","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"with:","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"the residuals r_imathbbR^nrightarrowmathbbR and the constraints c_imathbbR^nrightarrowmathbbR assumed to be mathcalC^2 functions;\nnorm cdot denoting the Euclidean norm.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Note that with this formulation, bounds constraints are not distinguished from general inequality constraints. Though, as shown later in this section, they can be provided as vectors of lower and/or upper bounds, which is more convenient for this type of constraints.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"It should be borne in mind, however, that the method implemented in Enlsip.jl has been conceived for nonlinear problems, as there is no other assumption made about the nature of the residuals and constraints functions, apart from being two-time continously differentiable. The algorithm can still be used to solve linear least squares subject to linear constraints but it will not be as effective as other software where those aspects are taken into account in the design of the optimization method.","category":"page"},{"location":"tutorial/#Instantiate-a-model","page":"Usage","title":"Instantiate a model","text":"","category":"section"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Solving a problem with Enlsip is organized in two steps.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"First, a model of type CnlsModel must be instantiated.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"The CnlsModel constructor requires the evaluation functions of residuals, constraints, their associated jacobian matrices and dimensions of the problem.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Although the package enables one to create linear unconstrained least squares, it is recommended to use it to solve nonlinear least squares with general constraints.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"The three following positional arguments are mandatory to create a model:","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"residuals : function that computes the vector of residuals\nnb_parameters : number of variables\nnb_residuals : number of residuals","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"The following keywords arguments are optional and deal with constraints and Jacobian matrices. If the Jacobian matrices functions are not provided, they are computed numerically by forward differences using automatic differenciation[Backend].","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"[Backend]: ForwardDiff.jl https://juliadiff.org/ForwardDiff.jl/stable/","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Argument Details\nstarting_point initial solution (can be an infeasbile point)\njacobian_residuals function computing the Jacobian matrix of the residuals\neq_constraints function computing the equality constraints\njacobian_eqcons function computing the Jacobian matrix of the equality constraints\nnb_eqcons number of equality constraints\nineq_constraints function computing the inequality constraints\njacobian_ineqcons function computing the Jacobian matrix of the inequality constraints\nnb_ineqcons number of inequality constraints\nx_low vector of lower bounds\nx_upp vector of upper bounds","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"It is assumed that the the different functions passed as arguments of the CnlsModel constructor are called as f(x), where x is a vector of nb_parameters elements and f is one of the functions residuals, eq_constraints, jacobian_eqcons etc.","category":"page"},{"location":"tutorial/#Solving-a-model","page":"Usage","title":"Solving a model","text":"","category":"section"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Then, the Enlsip solver can be used by calling the solve! function on a instantiated model. By default, the tolerance used in the algorithm is the square root of the relative precision on floating point numbers. For instance, with Float64, it will approximately equal 1e-8.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Enlsip.solve!","category":"page"},{"location":"tutorial/#Enlsip.solve!","page":"Usage","title":"Enlsip.solve!","text":"solve!(model)\n\nOnce a CnlsModel has been instantiated, this function solves the optimzation problem associated by using the method implemented in the Enlsip solver.\n\nKeywords arguments:\n\nsilent::Bool \nSet to false if one wants the algorithm to print details about the iterations and termination of the solver\nDefault value is true, i.e. by default, there is no output. If one wants to print those information afert solving, the print_cnls_model method \ncan be called.\nmax_iter::Int \nMaximum number of iterations allowed\nDefault value is set to 100\nscaling::Bool \nSet to true if one wants the algorithm to work with a constraints jacobian matrix whose rows are scaled (i.e. all constraints gradients vectors are scaled)\nDefault value is set to false\n\n\n\n\n\n","category":"function"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Diagnosis of the conduct of the algorithm can be printed by either setting the silent keyword argument of the function solve! to false or by calling print_cnls_model after solving. Here are some details on how to read and understand the different columns of the output:","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Column Description\niter iteration number\nobjective value of the sum of squared residuals (i.e. objective function) at current point\nvertvert active_constraints vertvert^2 value of the sum of squared active constraints at current point\nvertvert p vertvert norm of the search direction computed at current iteration\nalpha value of the steplength computed at current iteration\nreduction reduction in the objective function performed after moving to the next iterate","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"One can get additional info about termination of the algorithm by calling one of the following functions:","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Name\nsolution\nstatus\nconstraints_values\nobjective_value","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Enlsip.solution","category":"page"},{"location":"tutorial/#Enlsip.solution","page":"Usage","title":"Enlsip.solution","text":"solution(model)\n\nOnce the given model has been solved, this function returns the optimal solution, or last solution obtained if no convergence, as a Vector of approriate dimension.\n\n\n\n\n\n","category":"function"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Enlsip.status","category":"page"},{"location":"tutorial/#Enlsip.status","page":"Usage","title":"Enlsip.status","text":"status(model)\n\nThis functions returns a Symbol that gives brief information on the solving status of model.\n\nIf a model has been instantiated but the solver has not been called yet, it will return :unsolved.\n\nOnce the solver has been called and if a first order critical point satisfying the convergence criteria has been computed, it will return :successfully_solved.\n\nIf the algorithm met an abnormall termination criteria, it will return one of the following:\n\n:failed : the algorithm encoutered a numerical error that triggered termination\n:maximum_iterations_exceeded : a solution could not be reached within the maximum number of iterations.\n\n\n\n\n\n","category":"function"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Enlsip.constraints_values","category":"page"},{"location":"tutorial/#Enlsip.constraints_values","page":"Usage","title":"Enlsip.constraints_values","text":"constraints_values(model)\n\nComputes values of all the constraints in model at the solution. \n\nThe vector returned is the concatenation of equalities, inequalities and box constraints (in that order).\n\nFor instance, let xₛ be the solution found. If functions h, g compute equality and inequality constraints and xₗ, xᵤ are vectors of lower and lower bounds, it will return [h(xₛ); g(xₛ); xₛ-xₗ; xᵤ-xₛ].\n\nIf one wants to compute each type of constraints seperately, see equality_constraints_values, inequality_constraints_values and bounds_constraints_values.\n\n\n\n\n\n","category":"function"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Enlsip.objective_value","category":"page"},{"location":"tutorial/#Enlsip.objective_value","page":"Usage","title":"Enlsip.objective_value","text":"objective_value(model)\n\nOnce the given model has been solved, returns the value of the objective function, i.e. sum of squared residuals functions, computed at the optimal solution. If no convergence, this value is computed at the last solution obtained.\n\n\n\n\n\n","category":"function"},{"location":"tutorial/#Examples","page":"Usage","title":"Examples","text":"","category":"section"},{"location":"tutorial/#Problem-65-from-Hock-and-Schittkowski-collection[HS80]","page":"Usage","title":"Problem 65 from Hock and Schittkowski collection[HS80]","text":"","category":"section"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"We show how to implement and solve the following problem:","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"beginaligned\nmin_x_1 x_2 x_3 quad  (x_1-x_2)^2 + dfrac(x_1+x_2-10)^29+(x_3-5)^2  \ntextst quad  48-x_1^2-x_2^2-x_3^2 geq 0\n -45leq x_i leq 45 quad i=12\n -5 leq x_3  leq 5\nendaligned","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"The expected optimal solution is (3650461821 365046168 46204170507).","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Associated value of objective function equals 09535288567.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"First, we provide the dimensions of the problems.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"# Dimensions of the problem\n\nn = 3 # number of parameters\nm = 3 # number of residuals\nnb_eq = 0 # number of equality constraints\nnb_constraints = 7 # number of inequality constraints\nnothing # hide","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Then, we define the functions required to compute the residuals, constraints, their respective jacobian matrices and a starting point. In this example, we use the starting point given in the reference[HS80], i.e. (-5 5 0)","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"# Residuals and Jacobian matrix associated\nr(x::Vector) = [x[1] - x[2]; (x[1]+x[2]-10.0) / 3.0; x[3]-5.0]\n\njac_r(x::Vector) = [1. -1. 0;\n    1/3 1/3 0.;\n    0. 0. 1.]\n\n# Constraints (one equality and box constraints)\n\nc(x::Vector) = [48.0 - x[1]^2-x[2]^2-x[3]^2] # evaluation function for the equality constraint\njac_c(x::Vector) = [ -2x[1] -2x[2] -2x[3]] # Jacobian matrix of the equality constraint\n\nx_l = [-4.5, -4.5, -5.0] # lower bounds\nx_u = [4.5, 4.5, 5.0] # upper bounds\n\n# Starting point \nx0 = [-5.0, 5.0, 0.0]\nnothing # hide","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"A CnlsModel can now be instantiated.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"# Instantiate a model associated with the problem \nhs65_model = Enlsip.CnlsModel(r, n, m ;starting_point=x0, ineq_constraints = c, \nnb_ineqcons = 1, x_low=x_l, x_upp=x_u, jacobian_residuals=jac_r, jacobian_ineqcons=jac_c)\nnothing # hide","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Finally, the solve! function can be called on our model. In this example, keyword arguments remain to default values.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Enlsip.solve!(hs65_model)","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Once Enlsip solver has been executed on our problem, a summary of the conduct of the algorithm can be printed by calling print_cnls_model.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Enlsip.print_cnls_model(hs65_model)","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"If one just wants to know about termination of the algorithm, calling status will tell if the problem has been successfully solved or not.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Enlsip.status(hs65_model)","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"Then, calling solution and objective_value will respectively return the optimal solution obtained and the value of objective function at that point.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"hs65_solution = Enlsip.solution(hs65_model)","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"hs65_objective = Enlsip.objective_value(hs65_model)","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"The solution obtained is relatively close to the expected optimal solution, although it differs from more than the tolerance used.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"maximum(abs.(hs65_solution - [3.650461821, 3.65046168, 4.6204170507])) < sqrt(eps(Float64))","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"However, the difference between the objective value obtained with Enlsip and the expected one does not exceed the default tolerance.","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"abs(hs65_objective - 0.9535288567) < sqrt(eps(Float64))","category":"page"},{"location":"tutorial/","page":"Usage","title":"Usage","text":"[HS80]: W. Hock and K. Schittkowski. Test Examples for Nonlinear Programming Codes, volume 187 of Lecture Notes in Economics and Mathematical Systems. Springer, second edition, 1980.","category":"page"}]
}
