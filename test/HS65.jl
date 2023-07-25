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

    hs65_model = Enlsip.EnlsipModel(r, n, m ;jacobian_residuals=jac_r, starting_point=x0, ineq_constraints = c, jacobian_ineqcons=jac_c, nb_ineqcons = 1, x_low=x_l, x_upp=x_u)
    hs65_sol = solve(hs65_model,silent=true)

    @test r(x0) ≈ hs65_model.residuals.reseval(x0)
    @test vcat(c(x0),x0-x_l,x_u-x0) ≈ hs65_model.constraints.conseval(x0)
    @test jac_r(x0) ≈ hs65_model.residuals.jacres_eval(x0)
    @test nb_eq == hs65_model.nb_eqcons
    @test nb_constraints == hs65_model.nb_cons
    @test size(x0,1) == hs65_model.nb_parameters
    @test typeof(hs65_sol.exit_code) <: Int
    @test isfinite(hs65_sol.obj_value)
    @test size(hs65_sol.sol,1) == n
end