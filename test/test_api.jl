using Test
using Dedalus
using LinearAlgebra: norm

@testset "Package" begin
    @testset "Module loading" begin
        @test Dedalus.VERSION == "3.0.5"
    end

    @testset "Core types exist" begin
        @test isdefined(Dedalus, :Field)
        @test isdefined(Dedalus, :Distributor)
        @test isdefined(Dedalus, :CartesianCoordinates)
        @test isdefined(Dedalus, :Differentiate)
        @test isdefined(Dedalus, :Integrate)
        @test isdefined(Dedalus, :IVP)
        @test isdefined(Dedalus, :LBVP)
        @test isdefined(Dedalus, :CFL)
    end

    @testset "Milestone 2 types exist" begin
        @test isdefined(Dedalus, :MultidimensionalBasis)
        @test isdefined(Dedalus, :SpinBasis)
        @test isdefined(Dedalus, :PolarBasis)
        @test isdefined(Dedalus, :DiskBasis)
        @test isdefined(Dedalus, :AnnulusBasis)
        @test isdefined(Dedalus, :SphereBasis)
        @test isdefined(Dedalus, :PolarGradient)
        @test isdefined(Dedalus, :PolarDivergence)
        @test isdefined(Dedalus, :PolarLaplacian)
        @test isdefined(Dedalus, :MulCosine)
        @test isdefined(Dedalus, :SpinSkew)
        @test isdefined(Dedalus, :PolarTrace)
        @test isdefined(Dedalus, :RadialComponent)
        @test isdefined(Dedalus, :AngularComponent)
        @test isdefined(Dedalus, :AzimuthalComponent)
        @test isdefined(Dedalus, :SphereEllProduct)
        @test isdefined(Dedalus, :SWSHColatitudeTransform)
        @test isdefined(Dedalus, :DiskRadialTransform)
        @test isdefined(Dedalus, :SphereWrapper)
        @test isdefined(Dedalus, :Intertwiner)
    end

    @testset "Milestone 3 types exist" begin
        @test isdefined(Dedalus, :BallBasis)
        @test isdefined(Dedalus, :ShellBasis)
        @test isdefined(Dedalus, :SphericalCoordinates)
        @test isdefined(Dedalus, :SphericalGradient)
        @test isdefined(Dedalus, :SphericalDivergence)
        @test isdefined(Dedalus, :SphericalCurl)
        @test isdefined(Dedalus, :SphericalLaplacian)
        @test isdefined(Dedalus, :SphericalEllProduct)
        @test isdefined(Dedalus, :CrossProduct)
        @test isdefined(Dedalus, :DotProduct)
        @test isdefined(Dedalus, :CartesianTransposeComponents)
        @test isdefined(Dedalus, :CartesianTrace)
        @test isdefined(Dedalus, :RadialComponent)
        @test isdefined(Dedalus, :AngularComponent)
        @test isdefined(Dedalus, :S2_basis)
        @test hasfield(SphericalCoordinates, :S2coordsys)
    end
end

struct CopyNonSeparableTransform <: Dedalus.NonSeparableTransform end

function Dedalus.forward_reduced!(::CopyNonSeparableTransform, gdata, cdata)
    @test size(gdata) == size(cdata) == (1, 2, 5, 1)
    return copyto!(cdata, gdata)
end

function Dedalus.backward_reduced!(::CopyNonSeparableTransform, cdata, gdata)
    @test size(cdata) == size(gdata) == (1, 2, 5, 1)
    return copyto!(gdata, cdata)
end

@testset "Nonseparable transform coupled axes" begin
    plan = CopyNonSeparableTransform()
    input = reshape(collect(1:10), 2, 5)
    coeffs = zeros(Int, 2, 5)
    output = similar(coeffs)
    @test Dedalus.forward!(plan, input, coeffs, 1) == input
    @test Dedalus.backward!(plan, coeffs, output, 1) == input
end

@testset "Sphere degree maps" begin
    c = S2Coordinates("phi", "theta")
    d = Distributor(c, ComplexF64)
    b = SphereBasis(c, (8, 10), ComplexF64)
    maps = ell_maps(b, d)
    @test sort([entry[1] for entry in maps]) == collect(0:9)
    @test only(filter(entry -> entry[1] == 0, maps)) === (0, 1, 1)
    _, rows, columns = only(filter(entry -> entry[1] == 9, maps))
    domain = Dedalus.make_domain(d, (b,))
    m, ell = Dedalus._local_sphere_groups(b, Dedalus.coeff_layout(d), domain, (1, 1), 1)
    @test 9 in ell[rows, columns]

    coverage = zeros(Int, size(ell))
    for (degree, rows, columns) in Dedalus._ell_blocks(b, d)
        @test all(==(degree), ell[rows, columns])
        coverage[rows, columns] .+= 1
    end
    @test coverage == (@. abs(m) <= ell <= b.Lmax)
end

@testset "Colatitude inverse initializes omitted modes" for T in (Float64, ComplexF64)
    c = S2Coordinates("phi", "theta")
    d = Distributor(c, T)
    b = SphereBasis(c, (12, 5), T)
    plan = Dedalus.transform_plan(b, d, 5, 0)
    f = Field(d; bases = (b,))
    coeffs = f["c"]
    fill!(coeffs, zero(T))
    destination = ones(T, 12, 5)
    Dedalus.backward_reduced!(plan, reshape(coeffs, 1, size(coeffs)..., 1), reshape(destination, 1, 12, 5, 1))
    @test all(iszero, destination)
end

@testset "Colatitude inverse supports shared buffers" begin
    plan = SWSHColatitudeTransform(1, 0, ((0, 1:1, 1:1, 1:1),), 0)
    coeffs = fill(2.0, 1, 1)
    destination = similar(coeffs)
    backward!(plan, coeffs, destination, 1)
    @test destination[1] ≈ sqrt(2.0)
    backward!(plan, coeffs, coeffs, 1)
    @test coeffs[1] ≈ sqrt(2.0)
end

@testset "Ball radial transform supports integer selectors" begin
    c = SphericalCoordinates("phi", "theta", "r")
    d = Distributor(c, ComplexF64)
    b = BallBasis(c, (1, 1, 4), ComplexF64)
    maps = ell_maps(radial_basis(b), d)
    @test maps == ((0, 1, 1),)
    plan = BallRadialTransform((1, 1, 4), 4, 3, maps, (), 0, 0, 0.0)
    samples = reshape(ComplexF64[1, 2, 3, 4], 1, 1, 4)
    coeffs = similar(samples)
    destination = similar(samples)
    forward!(plan, samples, coeffs, 3)
    backward!(plan, coeffs, destination, 3)
    @test norm(destination - samples) <= 1.0e-12
end
