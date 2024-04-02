# Unitary tests for non exported functions in Enlsip.jl

@testset "Working set structure tests" begin
    # Test with 5 equality constraints and 5 inequality constraints
    q1, l1 = 5, 10
    wrkset1 = Enlsip.WorkingSet(q1, l1)
    
    # Initialization of the working  set
    @test length(wrkset1.active) == l1 && length(wrkset1.inactive) == l1 - q1
    @test wrkset1.t == q1

    a1, d1 = 7, 10
    Enlsip.add_constraint!(wrkset1, a1 - wrkset1.t)
    Enlsip.add_constraint!(wrkset1, d1- wrkset1.t)
    Enlsip.remove_constraint!(wrkset1, 7)

    @test a1 ∈ wrkset1.active && a1 ∉ wrkset1.inactive
    @test d1 ∉ wrkset1.active && d1 ∈ wrkset1.inactive
    @test wrkset1.t == q1 + 1
    @test count(!=(0), wrkset1.active) + count(!=(0), wrkset1.inactive) == l1

    # Test with no equality constraints and 8 inequality constraints
    q2, l2 = 0, 8
    wrkset2 = Enlsip.WorkingSet(q2,l2)
    
    @test wrkset2.active == zeros(l2) && wrkset2.inactive == [i for i ∈ 1:l2]
    # Make constraints 1, 4, 5 and 8 active
    Enlsip.add_constraint!(wrkset2, 1)
    Enlsip.add_constraint!(wrkset2, 4-wrkset2.t)
    Enlsip.add_constraint!(wrkset2, 5-wrkset2.t)
    Enlsip.add_constraint!(wrkset2, 8-wrkset2.t)
    lmt = wrkset2.l - wrkset2.t
    @test wrkset2.t == 4 && lmt == 4
    @test all(∈(wrkset2.active[1:wrkset2.t]), [1,4,5,8])
    @test all(∈(wrkset2.inactive[1:lmt]), [2,3,6,7])
end