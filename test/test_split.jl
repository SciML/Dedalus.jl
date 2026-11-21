using Test
using Dedalus

@testset "String and expression splitting" begin
    @test natural_sort(["x10", "x2", "x1"]) == ["x1", "x2", "x10"]
    @test split_call("f(x, y)") == ("f", ("x", "y"))

    coords = CartesianCoordinates("x")
    dist = Distributor(coords, Float64)
    u = Field(dist; name = "u")
    v = Field(dist; name = "v")
    dt = time_derivative(u)
    @test Dedalus.split(dt + v, Dedalus.TimeDerivative) === (dt, v)

    product = 2 * v
    containing, independent = Dedalus.split(product, Dedalus.TimeDerivative)
    @test containing === 0
    @test independent.args == product.args
end
