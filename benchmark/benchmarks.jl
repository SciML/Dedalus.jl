using Dedalus, BenchmarkTools

const SUITE = BenchmarkGroup()

# NOTE: Field(dist) currently fails on master — field.jl pushes the field into
# dist.fields::WeakKeyDict without a Pair, so the transform/evaluate surface is
# unreachable. The suite benchmarks the coordinate/distributor/basis/grid
# construction path that does work.

# =============================================================================
# Coordinate / distributor / basis construction
# =============================================================================

c = Coordinate("x")

SUITE["construct"] = BenchmarkGroup()

SUITE["construct"]["coordinate"] = @benchmarkable Coordinate("x")
SUITE["construct"]["distributor"] = @benchmarkable Distributor(
    $c, Float64; mesh = Int[]
)
SUITE["construct"]["chebyshev"] = @benchmarkable ChebyshevT(
    $c, 256, (0.0, 1.0)
)
SUITE["construct"]["fourier"] = @benchmarkable Fourier(
    $c, 128, (0.0, 2π); dtype = Float64
)

dist = Distributor(c, Float64; mesh = Int[])
basis = ChebyshevT(c, 256, (0.0, 1.0))

# =============================================================================
# Grid evaluation
# =============================================================================

SUITE["grid"] = BenchmarkGroup()

SUITE["grid"]["local_grid"] = @benchmarkable local_grid($basis, $dist, 1)
SUITE["grid"]["local_grids"] = @benchmarkable local_grids($basis, $dist, (1,))
