@testset "Problem Chained Rosenbrock" begin

    n = 1000
    m = 2(n-1)
    nb_eq = n-2
    nb_constraints = n-2

    function r(x::Vector)
        n = length(x)
        m = 2(n-1)
        rx = Vector(undef,m)
        rx[1:n-1] = [10(x[i]^2 - x[i+1]) for i=1:n-1]
        rx[n:m] = [x[k-n+1] - 1 for k=n:m]
        return rx
    end
    
    # Constraints

    function c(x::Vector)
        n = length(x)
        cx = [3x[k+1]^3 + 2x[k+2] - 5 + sin(x[k+1]-x[k+2])*sin(x[k+1]+x[k+2]) + 4x[k+1] - 
            x[k]*exp(x[k]-x[k+1]) - 3 for k=1:n-2]
        return cx
    end

    function jac_res(x::Vector)
        n = size(x,1)
        m = 2(n-1)
        J = zeros(eltype(x), (m,n))

        for i=1:n-1
            J[i,i] = 20x[i]
            J[i,i+1] = -10
        end

        for i=n:m
            J[i,i-n+1] = 1
        end
        return J
    end

    function jac_cons(x::Vector)
        n = size(x,1)
        A = zeros(eltype(x), (n-2,n))
        for k=1:n-2
            A[k,k] = -(x[k]+1) * exp(x[k]-x[k+1])
            A[k,k+1] = 9x[k+1]^2 + cos(x[k+1]-x[k+2])*sin(x[k+1]+x[k+2]) + sin(x[k+1]-x[k+2])*cos(x[k+1]+x[k+2]) + 4 + x[k]*exp(x[k]-x[k+1])
            A[k,k+2] = 2 - cos(x[k+1]-x[k+2])*sin(x[k+1]+x[k+2]) + sin(x[k+1]-x[k+2])*cos(x[k+1]+x[k+2])
        end
        return A
    end
    
    x0 = [(mod(i,2) == 1 ? -1.2 : 1.0) for i=1:n]

    CRmodel = CnlsModel(r,n,m ;starting_point=x0, jacobian_residuals=jac_res, eq_constraints=c, jacobian_eqcons=jac_cons, nb_eqcons=nb_eq)

    CR_sol = solve(CRmodel,silent=false)

    @test r(x0) ≈ CRmodel.residuals.reseval(x0)
    @test c(x0) ≈ CRmodel.constraints.conseval(x0)
    @test jac_res(x0) ≈ CRmodel.residuals.jacres_eval(x0)
    @test jac_cons(x0) ≈ CRmodel.constraints.jaccons_eval(x0)
    @test nb_eq == CRmodel.nb_eqcons
    @test nb_constraints == CRmodel.nb_cons
    @test size(x0,1) == CRmodel.nb_parameters
    @test typeof(CR_sol.solved) == Bool
    @test isfinite(CR_sol.obj_value)
    @test size(CR_sol.sol,1) == n
    
end