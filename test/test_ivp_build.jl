using Test
using Dedalus

@testset "IVP construction and field arithmetic" begin
    c = CartesianCoordinates("x", "y", "z")
    d = Distributor(c, Float64)
    xb = RealFourier(c.coords[1], 8, (0, 2pi))
    yb = RealFourier(c.coords[2], 8, (0, 2pi))
    zb = RealFourier(c.coords[3], 8, (0, 2pi))
    x, y, z = local_grids(d, xb, yb, zb; scales = 1)

    @testset "field multiplication" begin
        f = Field(d, bases = (xb, yb, zb))
        g = Field(d, bases = (xb, yb, zb))
        f["g"] = @. sin(x) * sin(y) * sin(z)
        g["g"] .= 2
        h = evaluate(f * g)
        @test isapprox(h["g"], 2 .* f["g"]; atol = 1.0e-12)
    end

    @testset "string IVP and solver build" begin
        u = Field(d; name = "u", bases = (xb, yb, zb))
        v = VectorField(d, c; name = "v", bases = (xb, yb, zb))
        u["g"] = @. sin(x) * sin(y) * sin(z)
        problem = IVP([u, v])
        add_equation!(problem, "dt(u) + u = 0")
        add_equation!(problem, "dt(v) + v = 0")
        @test length(get_equations(problem)) == 2
        solver = build_solver(problem, "RK222")
        @test solver !== nothing
    end

    @testset "interpolate call syntax" begin
        u = Field(d; name = "u", bases = (xb, yb, zb))
        u["g"] = @. sin(x) * sin(y) * sin(z)
        slice = evaluate(u(x = 0))
        @test size(slice["g"], 1) == 1
        @test all(isapprox.(slice["g"], 0; atol = 1.0e-12))
    end

    @testset "S2 and radial basis multiplication" begin
        sc = SphericalCoordinates("phi", "theta", "r")
        shell = ShellBasis(sc, (8, 10, 8), ComplexF64; radii = (1, 2))
        combined = Dedalus.basis_mul(S2_basis(shell), radial_basis(shell))
        @test combined isa ShellBasis
        @test combined.shape == (8, 10, 8)
        ball = BallBasis(sc, (8, 10, 8), ComplexF64; radius = 1.5)
        combined_ball = Dedalus.basis_mul(radial_basis(ball), S2_basis(ball))
        @test combined_ball isa BallBasis
        @test combined_ball.shape == (8, 10, 8)
        @test Dedalus.basis_mul(nothing, xb) === xb
        @test Dedalus.basis_mul(xb, nothing) === xb
    end
end
