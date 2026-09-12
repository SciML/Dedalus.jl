using Test
using Dedalus

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
