module SHGO


using StaticArrays
using NonlinearOptimizationTestFunctions
using SciMLBase
using LinearAlgebra
import NonlinearOptimizationTestFunctions: lb, ub

# Alias
const NOTF = NonlinearOptimizationTestFunctions

# Teile laden
include("types.jl")
include("cache.jl")
include("triangulation/kuhn.jl")

export analyze, SHGOResult

function analyze(tf::NOTF.TestFunction; kwargs...)
    x_start = NOTF.start(tf)
    f_val   = tf.f(x_start)
    
    # Wir bauen eine einfache Struktur, die mit 'Any' in SHGOResult kompatibel ist
    # Das vermeidet den MethodError von SciMLBase komplett
    dummy_sol = (u = x_start, objective = f_val, retcode = ReturnCode.Success)
    
    return SHGOResult(
        dummy_sol,      # global_minimum
        [dummy_sol],    # local_minima
        0               # num_basins (hier schlägt der Test gleich FEHL, nicht ERROR)
    )
end

end # module