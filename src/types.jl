using StaticArrays
using NonlinearOptimizationTestFunctions
using SciMLBase

struct Simplex{N}
    vertices::Vector{SVector{N, Float64}}
    indices::Vector{CartesianIndex{N}}
end

struct Region{N}
    simplices::Vector{Simplex{N}}
end

# Wir machen die Result-Struktur einfach und ohne N-Parameter
struct SHGOResult
    global_minimum::Any
    local_minima::Vector{Any}
    num_basins::Int
end

export Simplex, Region, SHGOResult