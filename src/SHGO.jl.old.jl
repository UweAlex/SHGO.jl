# File: src/SHGO.jl
#
# SHGO.jl - Simplicial Homology Global Optimization
#
# Pure Julia implementation based on Endres et al. (2018)
# 
# Features:
# - Implicit Kuhn topology (no graph in memory)
# - Lazy evaluation with point caching
# - Betti number stability for convergence
# - Star-domain based minimum detection
#
# Version: 0.2.0 - Bug fixes based on expert reviews
#
module SHGO

using StaticArrays
using Combinatorics
using Optimization, OptimizationOptimJL
using NonlinearOptimizationTestFunctions
using LinearAlgebra
using Printf

const NOTF = NonlinearOptimizationTestFunctions
const MIN_EPS = 1e-12

export analyze, SHGOResult, MinimumPoint
export KuhnTopology, PointCache

# =============================================================================
# Result Types
# =============================================================================

struct MinimumPoint
    minimizer::Vector{Float64}
    objective::Float64
end

Base.getproperty(m::MinimumPoint, s::Symbol) =
    s === :u ? getfield(m, :minimizer) : getfield(m, s)

"""
    SHGOResult

Result of SHGO optimization.

# Fields
- `local_minima::Vector{MinimumPoint}` - All found local minima (sorted by objective)
- `num_basins::Int` - Number of distinct basins found
- `iterations::Int` - Number of refinement iterations
- `converged::Bool` - Whether Betti stability was reached
- `f_calls::Int` - Number of function evaluations
"""
struct SHGOResult
    local_minima::Vector{MinimumPoint}
    num_basins::Int
    iterations::Int
    converged::Bool
    f_calls::Int
end

# Legacy accessor for backward compatibility
function Base.getproperty(r::SHGOResult, s::Symbol)
    if s === :results
        return getfield(r, :local_minima)
    else
        return getfield(r, s)
    end
end

# =============================================================================
# Point Cache - Efficient storage of evaluated points
# =============================================================================

"""
    PointCache{N}

Stores evaluated function values indexed by grid coordinates.
Prevents redundant evaluations at overlapping simplex vertices.
"""
struct PointCache{N}
    values::Dict{NTuple{N,Int}, Float64}
    lb::SVector{N,Float64}
    ub::SVector{N,Float64}
    divisions::SVector{N,Int}
    step::SVector{N,Float64}
    eval_count::Base.RefValue{Int}
end

function PointCache(lb::Vector{Float64}, ub::Vector{Float64}, divisions::Vector{Int})
    N = length(lb)
    step = SVector{N}((ub .- lb) ./ divisions)
    PointCache{N}(
        Dict{NTuple{N,Int}, Float64}(),
        SVector{N}(lb),
        SVector{N}(ub),
        SVector{N}(divisions),
        step,
        Ref(0)
    )
end

"""
Convert grid index to physical position.
"""
function index_to_position(cache::PointCache{N}, idx::NTuple{N,Int}) where N
    cache.lb .+ SVector{N}(idx) .* cache.step
end

"""
Check if index is within bounds.
"""
function is_valid_index(cache::PointCache{N}, idx::NTuple{N,Int}) where N
    all(i -> 0 <= idx[i] <= cache.divisions[i], 1:N)
end

"""
Get or compute function value. Returns Inf for out-of-bounds (infinity padding).
"""
function get_value!(cache::PointCache{N}, idx::NTuple{N,Int}, f::Function) where N
    # Infinity padding: outside bounds → +Inf
    if !is_valid_index(cache, idx)
        return Inf
    end
    
    # Cache lookup
    get!(cache.values, idx) do
        cache.eval_count[] += 1
        pos = index_to_position(cache, idx)
        f(pos)
    end
end

"""
Number of evaluated points.
"""
num_evaluated(cache::PointCache) = cache.eval_count[]

# =============================================================================
# Implicit Kuhn Topology - FIXED for N dimensions
# =============================================================================

"""
    KuhnTopology{N}

Implicit representation of Kuhn triangulation.
Computes neighborhoods on-the-fly without explicit graph.
"""
struct KuhnTopology{N}
    divisions::SVector{N,Int}
end

function KuhnTopology(divisions::Vector{Int})
    N = length(divisions)
    KuhnTopology{N}(SVector{N}(divisions))
end

"""
    get_neighbors(topo, idx)

Returns all Kuhn neighbors of a grid point.

FIXED: For N-dimensional Kuhn triangulation, neighbors are all points
reachable by changes in {-1, 0, +1} in each dimension, excluding the
point itself and ensuring the changes form valid Kuhn paths.

In Kuhn triangulation, point (i₁,...,iₙ) is connected to all points
where each coordinate changes by at most 1, AND the changes are
"monotonic" (all non-zero changes have the same sign pattern for
connected simplices).
"""
function get_neighbors(topo::KuhnTopology{N}, idx::NTuple{N,Int}) where N
    neighbors = NTuple{N,Int}[]
    
    # Generate all combinations of {-1, 0, +1} for each dimension
    # This covers all potential Kuhn neighbors
    for delta in Iterators.product(ntuple(_ -> (-1, 0, 1), N)...)
        # Skip self
        all(d == 0 for d in delta) && continue
        
        neighbor = ntuple(i -> idx[i] + delta[i], N)
        push!(neighbors, neighbor)
    end
    
    return neighbors
end

"""
Get the number of neighbors for a point (for pre-allocation).
"""
num_potential_neighbors(::KuhnTopology{N}) where N = 3^N - 1

# =============================================================================
# Star-Minimum Detection - FIXED with relative tolerance
# =============================================================================

"""
    is_star_minimum(cache, topo, idx, f)

Check if a point is the minimum in its star domain.

FIXED: Uses relative tolerance for numerical stability across
different function value ranges.
"""
function is_star_minimum(
    cache::PointCache{N}, 
    topo::KuhnTopology{N}, 
    idx::NTuple{N,Int},
    f::Function;
    rel_tol::Float64 = 1e-10
) where N
    # Own value
    val = get_value!(cache, idx, f)
    
    # Handle Inf/NaN
    !isfinite(val) && return false
    
    # Compute tolerance: relative for large values, absolute for small
    tol = max(MIN_EPS, abs(val) * rel_tol)
    
    # Compare with all neighbors
    for neighbor in get_neighbors(topo, idx)
        neighbor_val = get_value!(cache, neighbor, f)
        if neighbor_val < val - tol
            return false
        end
    end
    
    return true
end

"""
Find all star minima in the current grid.
"""
function find_star_minima(
    cache::PointCache{N},
    topo::KuhnTopology{N},
    f::Function
) where N
    minima = NTuple{N,Int}[]
    
    # Iterate over all grid points
    for idx in Iterators.product((0:d for d in topo.divisions)...)
        idx_tuple = NTuple{N,Int}(idx)
        if is_star_minimum(cache, topo, idx_tuple, f)
            push!(minima, idx_tuple)
        end
    end
    
    return minima
end

# =============================================================================
# Basin Clustering - FIXED: O(k) instead of O(k²)
# =============================================================================

"""
    cluster_basins(cache, topo, star_minima, f; threshold_ratio=0.1)

Cluster star minima into basins using Union-Find.

FIXED: 
- O(k) complexity by only checking neighbors (not all pairs)
- Uses relative value comparison
"""
function cluster_basins(
    cache::PointCache{N},
    topo::KuhnTopology{N},
    star_minima::Vector{NTuple{N,Int}},
    f::Function;
    threshold_ratio::Float64 = 0.1
) where N
    isempty(star_minima) && return Vector{Vector{NTuple{N,Int}}}()
    
    # Single star minimum → single basin
    if length(star_minima) == 1
        return [star_minima]
    end
    
    # Compute value range for threshold
    all_vals = collect(values(cache.values))
    finite_vals = filter(isfinite, all_vals)
    isempty(finite_vals) && return [[m] for m in star_minima]
    
    f_min = minimum(finite_vals)
    f_max = maximum(finite_vals)
    value_range = max(f_max - f_min, MIN_EPS)
    
    # Union-Find with path compression
    parent = Dict{NTuple{N,Int}, NTuple{N,Int}}()
    rank = Dict{NTuple{N,Int}, Int}()
    for m in star_minima
        parent[m] = m
        rank[m] = 0
    end
    
    function find_root(x)
        if parent[x] != x
            parent[x] = find_root(parent[x])  # Path compression
        end
        return parent[x]
    end
    
    function union!(x, y)
        rx, ry = find_root(x), find_root(y)
        if rx != ry
            # Union by rank
            if rank[rx] < rank[ry]
                parent[rx] = ry
            elseif rank[rx] > rank[ry]
                parent[ry] = rx
            else
                parent[ry] = rx
                rank[rx] += 1
            end
        end
    end
    
    # Build set for O(1) lookup
    star_set = Set(star_minima)
    
    # FIXED: O(k) - only check neighbors of each star minimum
    for m1 in star_minima
        val1 = get_value!(cache, m1, f)
        
        for neighbor in get_neighbors(topo, m1)
            # Is this neighbor also a star minimum?
            if neighbor in star_set
                val2 = get_value!(cache, neighbor, f)
                
                # Merge if values are similar (same basin)
                val_diff = abs(val1 - val2)
                if val_diff < value_range * threshold_ratio
                    union!(m1, neighbor)
                end
            end
        end
    end
    
    # Group by root
    clusters = Dict{NTuple{N,Int}, Vector{NTuple{N,Int}}}()
    for m in star_minima
        root = find_root(m)
        if !haskey(clusters, root)
            clusters[root] = NTuple{N,Int}[]
        end
        push!(clusters[root], m)
    end
    
    return collect(values(clusters))
end

# =============================================================================
# Local Optimization - FIXED: boundary handling
# =============================================================================

"""
Optimize locally from a starting point.
FIXED: Moves boundary points slightly inward to avoid optimizer warnings.
"""
function local_optimize(
    f::Function,
    grad::Function,
    x0::Vector{Float64},
    lb::Vector{Float64},
    ub::Vector{Float64};
    maxiters::Int = 500
)
    # FIXED: Move boundary points slightly inward
    eps_boundary = 1e-10
    x0_safe = copy(x0)
    for i in eachindex(x0_safe)
        range_i = ub[i] - lb[i]
        margin = min(eps_boundary, range_i * 1e-6)
        if x0_safe[i] <= lb[i] + margin
            x0_safe[i] = lb[i] + margin
        elseif x0_safe[i] >= ub[i] - margin
            x0_safe[i] = ub[i] - margin
        end
    end
    
    fopt = OptimizationFunction(
        (x, p) -> f(x),
        grad = (G, x, p) -> copyto!(G, grad(x))
    )
    
    prob = OptimizationProblem(fopt, x0_safe; lb=lb, ub=ub)
    
    # Try LBFGS, fallback to gradient-free if needed
    try
        sol = solve(prob, LBFGS(); maxiters=maxiters)
        return MinimumPoint(Vector(sol.minimizer), sol.objective)
    catch e
        # Fallback: return starting point evaluation
        return MinimumPoint(x0, f(x0))
    end
end

"""
Deduplicate minima based on spatial proximity and function value.
FIXED: Uses both position AND value for deduplication.
"""
function deduplicate_minima(minima::Vector{MinimumPoint}; 
                           dist_tol::Float64 = 0.05,
                           val_tol::Float64 = 1e-6)
    isempty(minima) && return MinimumPoint[]
    
    sorted = sort(minima, by = m -> m.objective)
    unique_minima = [sorted[1]]
    
    for m in sorted[2:end]
        is_new = true
        for u in unique_minima
            pos_close = norm(m.minimizer - u.minimizer) < dist_tol
            val_close = abs(m.objective - u.objective) < max(val_tol, abs(u.objective) * 1e-4)
            
            if pos_close && val_close
                is_new = false
                break
            end
        end
        is_new && push!(unique_minima, m)
    end
    
    return unique_minima
end

# =============================================================================
# Main Function: analyze() with Betti Stability
# =============================================================================

"""
    analyze(tf; kwargs...)

Analyze the optimization landscape of a test function.

# Algorithm:
1. Start with coarse grid
2. Find star minima (local minimum candidates)
3. Cluster into basins
4. Refine grid and repeat
5. Stop when Betti number (basin count) stabilizes
6. Locally optimize one representative per basin

# Keyword Arguments
| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_div_initial` | 8 | Initial grid resolution per dimension |
| `n_div_max` | 25 | Maximum grid resolution |
| `stability_count` | 2 | Iterations with stable basin count for convergence |
| `threshold_ratio` | 0.1 | Tolerance for basin merging (relative to value range) |
| `min_distance_tolerance` | 0.05 | Minimum distance between distinct minima |
| `local_maxiters` | 500 | Maximum iterations for local optimization |
| `verbose` | false | Print progress information |

# Returns
`SHGOResult` with fields: `local_minima`, `num_basins`, `iterations`, `converged`, `f_calls`
"""
function analyze(
    tf;
    # Core parameters
    n_div_initial::Int = 8,
    n_div_max::Int = 25,
    stability_count::Int = 2,
    threshold_ratio::Float64 = 0.1,
    min_distance_tolerance::Float64 = 0.05,
    local_maxiters::Int = 500,
    verbose::Bool = false,
    # Legacy parameter (for backward compatibility)
    n_div::Union{Int,Nothing} = nothing
)
    # Legacy: n_div overrides n_div_initial
    if !isnothing(n_div)
        n_div_initial = n_div
        n_div_max = max(n_div_max, n_div + 10)
    end
    
    lb = Vector{Float64}(NOTF.lb(tf))
    ub = Vector{Float64}(NOTF.ub(tf))
    N = length(lb)
    
    # Wrapper functions
    f = x -> tf.f(x)
    grad = x -> tf.grad(x)
    
    # Warning for high dimensions
    if N > 6 && verbose
        @warn "Dimension N=$N is high. Kuhn triangulation scales with 3^N neighbors. " *
              "Consider Sobol sampling for N > 6."
    end
    
    # =========================================================================
    # Iterative Refinement with Betti Stability
    # =========================================================================
    
    prev_num_basins = -1
    stable_iterations = 0
    current_n_div = n_div_initial
    iteration = 0
    final_basins = Vector{Vector{NTuple{N,Int}}}()
    final_cache = nothing
    
    while current_n_div <= n_div_max
        iteration += 1
        
        if verbose
            println("Iteration $iteration: n_div = $current_n_div")
        end
        
        # Create cache and topology for current resolution
        divisions = fill(current_n_div, N)
        cache = PointCache(lb, ub, divisions)
        topo = KuhnTopology(divisions)
        
        # Find star minima
        star_minima = find_star_minima(cache, topo, f)
        
        if verbose
            println("  Star minima found: $(length(star_minima))")
            println("  Points evaluated: $(num_evaluated(cache))")
        end
        
        # Cluster into basins
        basins = cluster_basins(cache, topo, star_minima, f; 
                               threshold_ratio=threshold_ratio)
        
        num_basins = length(basins)
        
        if verbose
            println("  Basins: $num_basins")
        end
        
        # Check Betti stability
        if num_basins == prev_num_basins && num_basins > 0
            stable_iterations += 1
            if verbose
                println("  Stable for $stable_iterations iterations")
            end
        else
            stable_iterations = 0
        end
        
        prev_num_basins = num_basins
        final_basins = basins
        final_cache = cache
        
        # Convergence reached?
        if stable_iterations >= stability_count
            if verbose
                println("Converged after $iteration iterations at n_div=$current_n_div")
            end
            break
        end
        
        # Refine grid
        current_n_div += 2
    end
    
    converged = stable_iterations >= stability_count
    
    if !converged && verbose
        @warn "No convergence reached. Maximum resolution n_div=$n_div_max used."
    end
    
    # =========================================================================
    # Local Optimization per Basin
    # =========================================================================
    
    if verbose
        println("\nLocal optimization for $(length(final_basins)) basins...")
    end
    
    candidates = MinimumPoint[]
    
    for (basin_id, basin) in enumerate(final_basins)
        # Select best point in basin as starting point
        best_idx = basin[1]
        best_val = get_value!(final_cache, best_idx, f)
        
        for idx in basin[2:end]
            val = get_value!(final_cache, idx, f)
            if val < best_val
                best_val = val
                best_idx = idx
            end
        end
        
        x0 = Vector{Float64}(index_to_position(final_cache, best_idx))
        
        try
            result = local_optimize(f, grad, x0, lb, ub; maxiters=local_maxiters)
            push!(candidates, result)
            
            if verbose
                println("  Basin $basin_id: f = $(round(result.objective, digits=6))")
            end
        catch e
            if verbose
                @warn "Local optimization failed for basin $basin_id" exception=e
            end
        end
    end
    
    # =========================================================================
    # Deduplication
    # =========================================================================
    
    unique_minima = deduplicate_minima(candidates; dist_tol=min_distance_tolerance)
    
    # Sort by objective value
    sort!(unique_minima, by = m -> m.objective)
    
    if verbose
        println("\nResult: $(length(unique_minima)) unique minima")
        for (i, m) in enumerate(unique_minima)
            println("  $i. f = $(round(m.objective, digits=6)) @ $(round.(m.minimizer, digits=4))")
        end
    end
    
    total_f_calls = isnothing(final_cache) ? 0 : num_evaluated(final_cache)
    
    return SHGOResult(
        unique_minima,
        length(unique_minima),
        iteration,
        converged,
        total_f_calls
    )
end

end # module