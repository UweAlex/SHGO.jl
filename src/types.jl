module SHGO

using StaticArrays
using NonlinearOptimizationTestFunctions  # Deine Bibliothek

# Kern-Typen
struct Simplex{N}
    vertices::NTuple{N+1, SVector{N, Float64}}  # Koordinaten
    indices::NTuple{N+1, CartesianIndex{N}}    # Gitter-Indizes für Cache
end

struct Region{N}
    simplices::Vector{Simplex{N}}
end

struct SHGOResult{N}
    global_minimum::OptimizationSolution
    local_minima::Vector{OptimizationSolution}
    num_basins::Int
    # Optional: graph, homology etc.
end

export solve, analyze, SHGOResult

end  # module