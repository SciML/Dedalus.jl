"""
Tests for parallel HDF5 output of various dimensionality.

Translated from dedalus/tests_parallel/test_output_parallel.py.
Requires MPI with at least 4 processes and mesh=(2,2).
"""

using Test
using MPI
using HDF5
using Dedalus

@testset "Virtual output across writes" begin
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    dir = MPI.bcast(rank == 0 ? mktempdir() : "", 0, comm)
    coords = CartesianCoordinates("x", "y", "z")
    dist = Distributor(coords, Float64; mesh = (2, 2))
    bases = Tuple(RealFourier(coord, 8, (0, 2pi)) for coord in coords.coords)
    field = Field(dist; bases, name = "u")
    @test evaluate(field) === field
    @test length(unique(string.((field(x = 0, z = 0), field(y = 0, z = 0))))) == 2

    ev = Dedalus.Evaluator(dist, Dict())
    output = add_file_handler(ev, joinpath(dir, "snapshots"); parallel = "virtual")
    add_task!(output, field; name = "u")
    x, y, z = local_grids(dist, bases...)
    for write in 1:2
        field["g"] = @. write * cos(x) * cos(y) * cos(z)
        evaluate_handlers(ev, [output])
    end

    MPI.Barrier(comm)
    h5open(joinpath(dir, "snapshots", "snapshots_s1.h5"), "r") do file
        sl = Dedalus.slices(dist.grid_layout, field.domain, field.scales)
        dset = file["tasks/u"]
        @test size(dset) == (2, 8, 8, 8)
        for write in 1:2
            @test dset[write, sl...] ≈ (@. write * cos(x) * cos(y) * cos(z)) atol = 1.0e-12
        end
    end
    sync() do
        if rank == 0
            rm(dir; recursive = true)
        end
    end
end

@testset "Output Parallel" begin

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    # Skip if fewer than 4 processes
    if nprocs < 4
        @warn "test_output_parallel requires at least 4 MPI processes; skipping."
        @test_broken false
        return
    end

    dtype_range = [Float64, ComplexF64]
    dealias_range = [1, 1.5]

    # Helper: build Fourier basis matching dtype
    function make_fourier(coord, T; size, bounds, dealias)
        if T == ComplexF64
            return ComplexFourier(coord, size, bounds; dealias = dealias)
        else
            return RealFourier(coord, size, bounds; dealias = dealias)
        end
    end

    # ========================================================================
    # test_cartesian_output
    # Default parallel file handler (virtual): per-process files.
    # Each rank reads its own file and compares to evaluated fields.
    # ========================================================================
    @testset "cartesian output T=$T dealias=$dealias output_scales=$output_scales" for
        T in dtype_range,
            dealias in dealias_range,
            output_scales in [0.5, 1, 1.5]
        Nx = Ny = Nz = 16
        Lx = Ly = Lz = 2 * pi
        # Bases
        c = CartesianCoordinates("x", "y", "z")
        d = Distributor(c, T; mesh = (2, 2))
        xb = make_fourier(c.coords[1], T; size = Nx, bounds = (0, Lx), dealias = dealias)
        yb = make_fourier(c.coords[2], T; size = Ny, bounds = (0, Ly), dealias = dealias)
        zb = make_fourier(c.coords[3], T; size = Nz, bounds = (0, Lz), dealias = dealias)
        x, y, z = local_grids(d, xb, yb, zb; scales = 1)
        # Fields
        u = Field(d, name = "u", bases = (xb, yb, zb), dtype = T)
        v = VectorField(d, c, name = "v", bases = (xb, yb, zb), dtype = T)
        u["g"] = @. sin(x) * sin(y) * sin(z)
        # Problem
        problem = IVP([u, v])
        add_equation!(problem, "dt(u) + u = 0")
        add_equation!(problem, "dt(v) + v = 0")
        solver = build_solver(problem, "RK222")
        # Output -- default parallel mode (per-process files)
        test_dir = "test_output"
        tasks = [
            u, u(x = 0), u(y = 0), u(z = 0),
            u(x = 0, y = 0), u(x = 0, z = 0), u(y = 0, z = 0),
            u(x = 0, y = 0, z = 0),
            v, v(x = 0), v(y = 0), v(z = 0),
            v(x = 0, y = 0), v(x = 0, z = 0), v(y = 0, z = 0),
            v(x = 0, y = 0, z = 0),
        ]
        output = add_file_handler(solver.evaluator, test_dir, iter = 1)
        for task in tasks
            add_task!(output, task, layout = "g", name = string(task), scales = output_scales)
        end
        evaluate_handlers(solver.evaluator, [output])
        # Check solution -- each rank reads its per-process file
        errors = Float64[]
        h5_path = joinpath(
            test_dir, "$(basename(test_dir))_s1",
            "$(basename(test_dir))_s1_p$(rank).h5"
        )
        @test isfile(h5_path)
        if isfile(h5_path)
            h5open(h5_path, "r") do file
                for task in tasks
                    task_eval = evaluate(task)
                    change_scales!(task_eval, output_scales)
                    data = task_eval["g"]
                    if !isempty(data)
                        dset = file["tasks"][string(task)]
                        task_saved = dset[1, ntuple(_ -> Colon(), ndims(data))...]
                        @test size(task_saved) == size(data)
                        push!(errors, maximum(abs, data .- task_saved))
                    end
                end
            end
        end
        # Cleanup
        sync() do
            if mpi_rank() == 0
                rm(test_dir, recursive = true, force = true)
            end
        end
        @test all(e -> e < 1.0e-12, errors)
    end

    # ========================================================================
    # test_cartesian_output_virtual
    # Virtual file handler with virtual_file=True equivalent.
    # Each rank reads from the joint virtual file, extracting local slices.
    # ========================================================================
    @testset "cartesian output virtual T=$T dealias=$dealias output_scales=$output_scales" for
        T in dtype_range,
            dealias in dealias_range,
            output_scales in [0.5]
        Nx = Ny = Nz = 16
        Lx = Ly = Lz = 2 * pi
        # Bases
        c = CartesianCoordinates("x", "y", "z")
        d = Distributor(c, T; mesh = (2, 2))
        xb = make_fourier(c.coords[1], T; size = Nx, bounds = (0, Lx), dealias = dealias)
        yb = make_fourier(c.coords[2], T; size = Ny, bounds = (0, Ly), dealias = dealias)
        zb = make_fourier(c.coords[3], T; size = Nz, bounds = (0, Lz), dealias = dealias)
        x, y, z = local_grids(d, xb, yb, zb; scales = 1)
        # Fields
        u = Field(d, name = "u", bases = (xb, yb, zb), dtype = T)
        v = VectorField(d, c, name = "v", bases = (xb, yb, zb), dtype = T)
        u["g"] = @. sin(x) * sin(y) * sin(z)
        # Problem
        problem = IVP([u, v])
        add_equation!(problem, "dt(u) + u = 0")
        add_equation!(problem, "dt(v) + v = 0")
        solver = build_solver(problem, "RK222")
        # Output -- virtual file mode
        test_dir = "test_output"
        tasks = [
            u, u(x = 0), u(y = 0), u(z = 0),
            u(x = 0, y = 0), u(x = 0, z = 0), u(y = 0, z = 0),
            u(x = 0, y = 0, z = 0),
            v, v(x = 0), v(y = 0), v(z = 0),
            v(x = 0, y = 0), v(x = 0, z = 0), v(y = 0, z = 0),
            v(x = 0, y = 0, z = 0),
        ]
        output = add_file_handler(
            solver.evaluator, test_dir, iter = 1,
            max_writes = 1, parallel = "virtual"
        )
        for task in tasks
            add_task!(output, task, layout = "g", name = string(task), scales = output_scales)
        end
        evaluate_handlers(solver.evaluator, [output])
        # Check solution -- read from joint virtual file
        errors = Float64[]
        MPI.Barrier(comm)
        joint_path = joinpath(test_dir, "$(basename(test_dir))_s1.h5")
        @test isfile(joint_path)
        if isfile(joint_path)
            h5open(joint_path, "r") do file
                for task in tasks
                    task_name = string(task)
                    task_eval = evaluate(task)
                    change_scales!(task_eval, output_scales)
                    # Compute local slices from the global dataset
                    gl = d.grid_layout
                    local_sl = Dedalus.slices(gl, task_eval.domain, task_eval.scales)
                    tensor_dims = ntuple(i -> Colon(), length(task_eval.tensorsig))
                    full_slices = (tensor_dims..., local_sl...)
                    task_local = file["tasks"][task_name][1, full_slices...]
                    local_error = task_eval["g"] .- task_local
                    if length(local_error) > 0
                        push!(errors, maximum(abs.(local_error)))
                    end
                end
            end
        end
        # Cleanup
        sync() do
            if mpi_rank() == 0
                rm(test_dir, recursive = true, force = true)
            end
        end
        @test all(e -> e < 1.0e-12, errors)
    end

    # ========================================================================
    # test_cartesian_output_merged_virtual
    # Virtual output with post.merge_virtual_analysis.
    # ========================================================================
    @testset "cartesian output merged virtual T=$T dealias=$dealias output_scales=$output_scales" for
        T in dtype_range,
            dealias in dealias_range,
            output_scales in [0.5, 1]
        Nx = Ny = Nz = 16
        Lx = Ly = Lz = 2 * pi
        # Bases
        c = CartesianCoordinates("x", "y", "z")
        d = Distributor(c, T; mesh = (2, 2))
        xb = make_fourier(c.coords[1], T; size = Nx, bounds = (0, Lx), dealias = dealias)
        yb = make_fourier(c.coords[2], T; size = Ny, bounds = (0, Ly), dealias = dealias)
        zb = make_fourier(c.coords[3], T; size = Nz, bounds = (0, Lz), dealias = dealias)
        x, y, z = local_grids(d, xb, yb, zb; scales = 1)
        # Fields
        u = Field(d, name = "u", bases = (xb, yb, zb), dtype = T)
        v = VectorField(d, c, name = "v", bases = (xb, yb, zb), dtype = T)
        u["g"] = @. sin(x) * sin(y) * sin(z)
        # Problem
        problem = IVP([u, v])
        add_equation!(problem, "dt(u) + u = 0")
        add_equation!(problem, "dt(v) + v = 0")
        solver = build_solver(problem, "RK222")
        # Output -- virtual file mode
        test_dir = "test_output"
        tasks = [
            u, u(x = 0), u(y = 0), u(z = 0),
            u(x = 0, y = 0), u(x = 0, z = 0), u(y = 0, z = 0),
            u(x = 0, y = 0, z = 0),
            v, v(x = 0), v(y = 0), v(z = 0),
            v(x = 0, y = 0), v(x = 0, z = 0), v(y = 0, z = 0),
            v(x = 0, y = 0, z = 0),
        ]
        output = add_file_handler(
            solver.evaluator, test_dir, iter = 1,
            max_writes = 1, parallel = "virtual"
        )
        for task in tasks
            add_task!(output, task, layout = "g", name = string(task), scales = output_scales)
        end
        evaluate_handlers(solver.evaluator, [output])
        # Merge virtual datasets into a single file
        # (Julia equivalent of post.merge_virtual_analysis)
        if isdefined(Dedalus, :merge_virtual_analysis)
            Dedalus.merge_virtual_analysis(test_dir; cleanup = true)
        end
        MPI.Barrier(comm)
        # Check solution -- read from merged joint file
        errors = Float64[]
        joint_path = joinpath(test_dir, "$(basename(test_dir))_s1.h5")
        @test isfile(joint_path)
        if isfile(joint_path)
            h5open(joint_path, "r") do file
                for task in tasks
                    task_name = string(task)
                    task_eval = evaluate(task)
                    change_scales!(task_eval, output_scales)
                    # Compute local slices from the global dataset
                    gl = d.grid_layout
                    local_sl = Dedalus.slices(gl, task_eval.domain, task_eval.scales)
                    tensor_dims = ntuple(i -> Colon(), length(task_eval.tensorsig))
                    full_slices = (tensor_dims..., local_sl...)
                    task_local = file["tasks"][task_name][1, full_slices...]
                    local_error = task_eval["g"] .- task_local
                    if length(local_error) > 0
                        push!(errors, maximum(abs.(local_error)))
                    end
                end
            end
        end
        # Cleanup
        sync() do
            if mpi_rank() == 0
                rm(test_dir, recursive = true, force = true)
            end
        end
        @test all(e -> e < 1.0e-12, errors)
    end

    # ========================================================================
    # test_cartesian_output_merged
    # Per-process output with post.merge_analysis.
    # ========================================================================
    @testset "cartesian output merged T=$T dealias=$dealias output_scales=$output_scales" for
        T in dtype_range,
            dealias in dealias_range,
            output_scales in [0.5, 1]
        Nx = Ny = Nz = 16
        Lx = Ly = Lz = 2 * pi
        # Bases
        c = CartesianCoordinates("x", "y", "z")
        d = Distributor(c, T; mesh = (2, 2))
        xb = make_fourier(c.coords[1], T; size = Nx, bounds = (0, Lx), dealias = dealias)
        yb = make_fourier(c.coords[2], T; size = Ny, bounds = (0, Ly), dealias = dealias)
        zb = make_fourier(c.coords[3], T; size = Nz, bounds = (0, Lz), dealias = dealias)
        x, y, z = local_grids(d, xb, yb, zb; scales = 1)
        # Fields
        u = Field(d, name = "u", bases = (xb, yb, zb), dtype = T)
        v = VectorField(d, c, name = "v", bases = (xb, yb, zb), dtype = T)
        u["g"] = @. sin(x) * sin(y) * sin(z)
        # Problem
        problem = IVP([u, v])
        add_equation!(problem, "dt(u) + u = 0")
        add_equation!(problem, "dt(v) + v = 0")
        solver = build_solver(problem, "RK222")
        # Output -- per-process file mode (not virtual)
        test_dir = "test_output"
        tasks = [
            u, u(x = 0), u(y = 0), u(z = 0),
            u(x = 0, y = 0), u(x = 0, z = 0), u(y = 0, z = 0),
            u(x = 0, y = 0, z = 0),
            v, v(x = 0), v(y = 0), v(z = 0),
            v(x = 0, y = 0), v(x = 0, z = 0), v(y = 0, z = 0),
            v(x = 0, y = 0, z = 0),
        ]
        output = add_file_handler(
            solver.evaluator, test_dir, iter = 1,
            max_writes = 1, parallel = "virtual"
        )
        for task in tasks
            add_task!(output, task, layout = "g", name = string(task), scales = output_scales)
        end
        evaluate_handlers(solver.evaluator, [output])
        # Merge per-process files into a single file
        # (Julia equivalent of post.merge_analysis)
        if isdefined(Dedalus, :merge_analysis)
            Dedalus.merge_analysis(test_dir; cleanup = true)
        end
        MPI.Barrier(comm)
        # Check solution -- read from merged file
        errors = Float64[]
        merged_path = joinpath(test_dir, "$(basename(test_dir))_s1.h5")
        @test isfile(merged_path)
        if isfile(merged_path)
            h5open(merged_path, "r") do file
                for task in tasks
                    task_name = string(task)
                    task_eval = evaluate(task)
                    change_scales!(task_eval, output_scales)
                    # Compute local slices from the global dataset
                    gl = d.grid_layout
                    local_sl = Dedalus.slices(gl, task_eval.domain, task_eval.scales)
                    tensor_dims = ntuple(i -> Colon(), length(task_eval.tensorsig))
                    full_slices = (tensor_dims..., local_sl...)
                    task_local = file["tasks"][task_name][1, full_slices...]
                    local_error = task_eval["g"] .- task_local
                    if length(local_error) > 0
                        push!(errors, maximum(abs.(local_error)))
                    end
                end
            end
        end
        # Cleanup
        sync() do
            if mpi_rank() == 0
                rm(test_dir, recursive = true, force = true)
            end
        end
        @test all(e -> e < 1.0e-12, errors)
    end

end
