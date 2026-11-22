using Test
using MPI
using Dedalus

@testset "Distributed layout transitions" begin
    @testset "Fourier transpose T=$T" for T in (Float64, ComplexF64)
        coords = CartesianCoordinates("x", "y", "z")
        dist = Distributor(coords, T; mesh = (2, 2))
        fourier = T === Float64 ? RealFourier : ComplexFourier
        bases = Tuple(fourier(coord, n, (0, 2pi)) for (coord, n) in zip(coords.coords, (10, 14, 18)))
        x, y, z = local_grids(dist, bases...)
        field = Field(dist; bases, dtype = T)
        values = @. T(sin(x) * cos(2y) * sin(3z))
        field["g"] = values

        Dedalus.change_layout!(field, "c")
        Dedalus.change_layout!(field, "g")
        @test field["g"] ≈ values atol = 1.0e-12
    end

    @testset "Evaluator traverses every layout" begin
        coords = CartesianCoordinates("x", "y", "z")
        dist = Distributor(coords, Float64; mesh = (2, 2))
        bases = Tuple(RealFourier(coord, 8, (0, 2pi)) for coord in coords.coords)
        field = Field(dist; bases)
        x, y, z = local_grids(dist, bases...)
        expected = @. cos(x) + 0 * y + 0 * z
        field["g"] = expected
        ev = Dedalus.Evaluator(dist, Dict())
        Dedalus.require_coeff_space(ev, [field])
        @test field.layout === dist.coeff_layout
        Dedalus.require_grid_space(ev, [field])
        @test field.layout === dist.grid_layout
        @test field["g"] ≈ expected atol = 1.0e-12
    end
end
