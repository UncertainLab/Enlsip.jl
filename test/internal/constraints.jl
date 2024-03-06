@testset "ConstraintsFunction structure" begin

    
    c(x::Vector) = [3x[1]^3 + 2x[2] - 5 + sin(x[1]-x[2]*sin(x[1]+x[2])), 
    4x[4] - x[3]*exp(x[3]-x[4]) - 3]
    x = zeros(4)
    c_func = Enlsip.ConstraintsFunction(c)
    @test c_func.nb_conseval == 0 && c_func.nb_jaccons_eval == 0
    A = Matrix{Float64}(undef, 2,4)
    Enlsip.jaccons_eval!(c_func, x,A)
    @test c_func.nb_conseval == 0 && c_func.nb_jaccons_eval == 1

    x_low = [-1.0, -Inf, -2.0, -Inf]
    x_upp = [Inf, Inf, 5.0, 10.0]
    f_bounds, jac_f_bounds = Enlsip.box_constraints(x_low, x_upp)
    @test Inf ∉ abs.(f_bounds(x))
    @test eltype(jac_f_bounds(x)) <: AbstractFloat && all(isfinite, jac_f_bounds(x))

    h_func = Enlsip.instantiate_constraints_w_bounds(c, nothing, nothing, nothing, x_low, x_upp)
    hx = Vector{Float64}(undef, 6)
    Ah = Matrix{Float64}(undef, 6, 4)
    Enlsip.cons_eval!(h_func, x, hx)
    Enlsip.jaccons_eval!(h_func, x, Ah)
    @test hx ≈ h_func.conseval(x)
    @test Ah ≈ h_func.jaccons_eval(x) 
end