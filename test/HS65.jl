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


    hs65_model = Enlsip.CnlsModel(r, n, m ;jacobian_residuals=jac_r, starting_point=x0, ineq_constraints = c, jacobian_ineqcons=jac_c, nb_ineqcons = 1, x_low=x_l, x_upp=x_u)
    solve!(hs65_model, silent=false)


    @test size(x0,1) == hs65_model.nb_parameters
    @test r(x0) ≈ hs65_model.residuals(x0)
    @test jac_r(x0) ≈ hs65_model.jacobian_residuals(x0)
    @test vcat(c(x0),x0-x_l,x_u-x0) ≈ vcat(hs65_model.ineq_constraints(x0), x0-hs65_model.x_low, hs65_model.x_upp-x0)
    @test nb_constraints == total_nb_constraints(hs65_model)
    @test status(hs65_model) in values(dict_status_codes)
    @test typeof(solution(hs65_model)) <: Vector && size(solution(hs65_model),1) == n
    @test typeof(objective_value(hs65_model)) <: Number && isfinite(objective_value(hs65_model))


end