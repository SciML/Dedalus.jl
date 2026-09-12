using Test
using Dedalus
using SciMLTesting
using JET: JET

@testset "Quality Assurance" begin
    run_qa(
        Dedalus;
        aqua_kwargs = (;
            # MPI.jl is a deliberate lazy dependency: it is loaded on demand
            # through `Base.require(Main, :MPI)` in `init_mpi!`, never imported
            # during module load.
            stale_deps = (; ignore = [:MPI]),
        ),
        ei_kwargs = (;
            all_qualified_accesses_are_public = (;
                ignore = (
                    # MPI is loaded lazily via Base.require(Main, :MPI) in init_mpi!
                    :require,
                    # HDF5.API internals are required for parallel/virtual datasets
                    :API, :Dataset, :H5P_DATASET_CREATE, :H5P_DEFAULT, :H5S_SELECT_SET,
                    :h5d_create, :h5l_move, :h5p_close, :h5p_create, :h5p_set_virtual,
                    :h5s_close, :h5s_create_simple, :h5s_select_hyperslab,
                    :has_parallel, :set_extent_dims,
                    # FFTW planner flags are not exported by FFTW.jl
                    :ESTIMATE, :EXHAUSTIVE, :MEASURE, :PATIENT,
                    # Banded LU factorisations have no public LAPACK spelling
                    :gbtrf!, :gbtrs!,
                    # DedalusLogger needs Logging's internal filtering interface
                    :catch_exceptions, :handle_message, :min_enabled_level, :shouldlog,
                    # Equation parsing caches Base method internals
                    :eval, :method_argnames,
                    # `Main.get_basis_axis` is a deliberate extension point for
                    # user code that predates the owner-module convention.
                    :get_basis_axis,
                ),
            ),
            # `Main.get_basis_axis` is a deliberate extension point for user
            # code that predates the owner-module convention.
            all_qualified_accesses_via_owners = (; ignore = (:get_basis_axis,)),
        ),
    )
end
