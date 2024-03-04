# Unitary tests for non exported functions in Enlsip.jl

@testset "Working set structure tests" begin
    # Test with 5 equality constraints and 5 inequality constraints
    q1, l1 = 5, 10
    test_wrkset1 = Enlsip.WorkingSet(q1, l1)
    
    # Initialization of the working  set
    @test length(test_wrkset1.active) == l1 && length(test_wrkset1.inactive) == l1 - q1
    @test test_wrkset1.t == q1

    a1, d1 = 7, 10
    Enlsip.add_constraint!(test_wrkset1, a1 - test_wrkset1.t)
    Enlsip.add_constraint!(test_wrkset1, d1- test_wrkset1.t)
    Enlsip.remove_constraint!(test_wrkset1, 7)

    @test a1 ∈ test_wrkset1.active && a1 ∉ test_wrkset1.inactive
    @test d1 ∉ test_wrkset1.active && d1 ∈ test_wrkset1.inactive
    @test test_wrkset1.t == q1 + 1
    @test count(!=(0), test_wrkset1.active) + count(!=(0), test_wrkset1.inactive) == l1

    # Test with no equality constraints and 8 inequality constraints
    q2, l2 = 0, 8
    test_wrkset2 = Enlsip.WorkingSet(q2,l2)
    
    @test test_wrkset2.active == zeros(l2) && test_wrkset2.inactive == [i for i ∈ 1:l2]
    # Make constraints 1, 4, 5 and 8 active
    Enlsip.add_constraint!(test_wrkset2, 1)
    Enlsip.add_constraint!(test_wrkset2, 4-test_wrkset2.t)
    Enlsip.add_constraint!(test_wrkset2, 5-test_wrkset2.t)
    Enlsip.add_constraint!(test_wrkset2, 8-test_wrkset2.t)

    @test test_wrkset2.t == 4
    @test test_wrkset2.active[1:test_wrkset2.t] == [1, 4, 5, 8]
    @test test_wrkset2.inactive[1:test_wrkset2.l - test_wrkset2.t] == [2, 3, 6, 7]
end