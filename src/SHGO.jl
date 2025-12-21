module SHGO

using StaticArrays
using LazySets
using Graphs
using Optimization
using OptimizationOptimJL
using ConcurrentCollections
using NonlinearOptimizationTestFunctions

include("types.jl")
include("cache.jl")
include("triangulation/kuhn.jl")
include("pruning/gradient_hull.jl")
include("pruning/value_pruning.jl")
include("clustering.jl")
include("local_search.jl")

export analyze, SHGOResult

# Haupt-API – später erweitern
function analyze(tf::TestFunction; kwargs...)
    # Platzhalter – wird in Phase 1 gefüllt
    error("Not implemented yet – Phase 1 in progress")
end

end # module SHGO