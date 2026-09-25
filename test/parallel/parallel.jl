using Test
using Dedalus
using MPI

@testset "Parallel (MPI)" begin
    mpi_dir = joinpath(@__DIR__, "mpi")
    project = Base.active_project()

    exe = MPI.mpiexec()
    for f in sort(readdir(mpi_dir))
        endswith(f, ".jl") || continue
        file = joinpath(mpi_dir, f)
        @testset "$f" begin
            script = tempname() * ".jl"
            write(
                script, """
                using MPI
                MPI.Init()
                using Dedalus
                Dedalus.init_mpi!()
                using Test
                include($(repr(file)))
                MPI.Finalize()
                """
            )
            cmd = `$exe -n 4 $(Base.julia_cmd()) --project=$project $script`
            @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
        end
    end
end
