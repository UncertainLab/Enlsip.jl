@testset "Osborne 2: Box constrained problem" begin
    
    n = 11
    m = 65 
    nb_eq = 0
    nb_constraints = 22

    # DataPoints

    dataset = [1 0.0 1.366;
    2 0.1 1.191;
    3 0.2 1.112;
    4 0.3 1.013;
    5 0.4 0.991;
    6 0.5 0.885;
    7 0.6 0.831;
    8 0.7 0.847;
    9 0.8 0.786;
    10 0.9 0.725;
    11 1.0 0.746;
    12 1.1 0.679;
    13 1.2 0.608;
    14 1.3 0.655;
    15 1.4 0.616;
    16 1.5 0.606;
    17 1.6 0.602;
    18 1.7 0.626;
    19 1.8 0.651;
    20 1.9 0.724;
    21 2.0 0.649;
    22 2.1 0.649;
    23 2.2 0.694;
    24 2.3 0.644;
    25 2.4 0.624;
    26 2.5 0.661;
    27 2.6 0.612;
    28 2.7 0.558;
    29 2.8 0.533;
    30 2.9 0.495;
    31 3.0 0.500;
    32 3.1 0.423;
    33 3.2 0.395;
    34 3.3 0.375;
    35 3.4 0.538;
    36 3.5 0.522;
    37 3.6 0.506;
    38 3.7 0.490;
    39 3.8 0.478;
    40 3.9 0.467;
    41 4.0 0.457;
    42 4.1 0.457;
    43 4.2 0.457;
    44 4.3 0.457;
    45 4.4 0.457;
    46 4.5 0.457;
    47 4.6 0.457;
    48 4.7 0.457;
    49 4.8 0.457;
    50 4.9 0.457;
    51 5.0 0.457;
    52 5.1 0.431;
    53 5.2 0.431;
    54 5.3 0.424;
    55 5.4 0.420;
    56 5.5 0.414;
    57 5.6 0.411;
    58 5.7 0.406;
    59 5.8 0.406;
    60 5.9 0.406;
    61 6.0 0.406;
    62 6.1 0.406;
    63 6.2 0.406;
    64 6.3 0.406;
    65 6.4 0.406]

    t = dataset[1:m,2]
    y = dataset[1:m,3]

    

    function r_k(x::Vector, t::T, y::T) where {T}
        rx = x[1]*exp(-x[5]*t) + x[2]*exp(-x[6]*(t-x[9])^2) + x[3]*exp(-x[7]*(t-x[10])^2) + x[4]*exp(-x[8]*(t-x[11])^2)
        return y - rx
    end
    
    r(x::Vector) = [r_k(x,t[k],y[k]) for k=1:m]

    low_bounds = [1.31, 0.4314, 0.6336, 0.5, 0.5, 0.6, 1.0, 4.0, 2.0, 4.5689, 5.0]
    upp_bounds = [1.4, 0.8, 1.0, 1.0, 1.0, 3.0, 5.0, 7.0, 2.5, 5.0, 6.0]

    # Saved starting point
    x0 = [1.3344098963722457
    0.5572842161127423
    0.6757364753061974
    0.8291980513226953
    0.9233565833014519
    0.9588470511477797
    1.9610314699563896
    4.055321823656234
    2.048625993866472
    4.60296578920499
    5.95212572157736]

    osborne2_model = Enlsip.CnlsModel(r,n,m; starting_point = x0, x_low = low_bounds, x_upp = upp_bounds)
    solve!(osborne2_model)
    bounds_values = Enlsip.constraints_values(osborne2_model)
    @test nb_upper_bounds(osborne2_model) == length(upp_bounds) && nb_lower_bounds(osborne2_model) == length(low_bounds)
    @test total_nb_constraints(osborne2_model) == nb_constraints
    @test osborne2_model.jacobian_residuals === nothing
    @test nb_equality_constraints(osborne2_model) == 0 && osborne2_model.eq_constraints === nothing
    @test nb_inequality_constraints(osborne2_model) == 0 && osborne2_model.ineq_constraints === nothing
    @test length(bounds_values) == length(low_bounds) + length(upp_bounds)
end