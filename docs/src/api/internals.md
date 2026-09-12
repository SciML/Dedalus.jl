# [Internals](@id Internals)

Low-level machinery: caching and dispatch infrastructure, sparse-matrix
utilities, Jacobi/Clenshaw kernels, MPI helpers, equation parsing, and the
coefficient/field system containers.

## Field and coefficient systems

```@autodocs
Modules = [Dedalus]
Pages = ["core/system.jl"]
```

## Caching and dispatch

```@autodocs
Modules = [Dedalus]
Pages = ["tools/cache.jl", "tools/dispatch.jl"]
```

## Exceptions

```@autodocs
Modules = [Dedalus]
Pages = ["tools/exceptions.jl"]
```

## Array and linear-algebra utilities

```@autodocs
Modules = [Dedalus]
Pages = ["tools/array.jl", "tools/linalg.jl", "tools/jacobi.jl", "tools/clenshaw.jl"]
```

## General utilities and parsing

```@autodocs
Modules = [Dedalus]
Pages = ["tools/general.jl", "tools/parsing.jl", "tools/logging.jl",
         "tools/random_arrays.jl", "tools/config.jl"]
```

## MPI and parallelism

```@autodocs
Modules = [Dedalus]
Pages = ["tools/parallel.jl"]
```
