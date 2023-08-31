@testset "Problem 65 from Hock Schittkowski collection" begin
    n = 3
    m = 3
    nb_eq = 0
    nb_constraints = 7

    r(x::Vector) = [x[1] - x[2]; (x[1]+x[2]-10.0) / 3.0; x[3]-5.0]
    jac_r(x::Vector) = [1. -1. 0;
    1/3 1/3 0.;
    0. 0. 1.]

    c(x::Vector) = [48.0 - x[1]^2-x[2]^2-x[3]^2]
    jac_c(x::Vector) = [ -2x[1] -2x[2] -2x[3]]
    x_l = [-4.5, -4.5, -5.0]
    x_u = [4.5, 4.5, 5.0] 

    x0 = [-5.0, 5.0, 0.0]

    hs65_model = Enlsip.CNLSModel(r, n, m ;jacobian_residuals=jac_r, starting_point=x0, ineq_constraints = c, jacobian_ineqcons=jac_c, nb_ineqcons = 1, x_low=x_l, x_upp=x_u)
    hs65_sol = solve(hs65_model;silent=true)

    hs65_model_2 = Enlsip.CnlsModel(r, n, m ;jacobian_residuals=jac_r, starting_point=x0, ineq_constraints = c, jacobian_ineqcons=jac_c, nb_ineqcons = 1, x_low=x_l, x_upp=x_u)
    solve!(hs65_model_2;silent=true)


    @test size(x0,1) == hs65_model_2.nb_parameters
    @test r(x0) ≈ hs65_model_2.residuals(x0)
    @test jac_r(x0) ≈ hs65_model_2.jacobian_residuals(x0)
    @test vcat(c(x0),x0-x_l,x_u-x0) ≈ vcat(hs65_model_2.ineq_constraints(x0), x0-hs65_model_2.x_low, hs65_model_2.x_upp-x0)
    @test nb_constraints == total_nb_constraints(hs65_model_2)
    @test status(hs65_model_2) in values(status_codes)
    @test typeof(solution(hs65_model_2)) <: Vector && size(solution(hs65_model_2),1) == n
    @test typeof(objective_value(hs65_model_2)) <: Number && isfinite(objective_value(hs65_model_2))

    @test r(x0) ≈ hs65_model.residuals.reseval(x0)
    @test vcat(c(x0),x0-x_l,x_u-x0) ≈ hs65_model.constraints.conseval(x0)
    @test jac_r(x0) ≈ hs65_model.residuals.jacres_eval(x0)
    @test nb_eq == hs65_model.nb_eqcons
    @test nb_constraints == hs65_model.nb_cons
    @test size(x0,1) == hs65_model.nb_parameters
    @test typeof(hs65_sol.solved) == Bool
    @test isfinite(hs65_sol.obj_value)
    @test size(hs65_sol.sol,1) == n
end