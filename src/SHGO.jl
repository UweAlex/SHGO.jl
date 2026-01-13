module SHGO

using StaticArrays
using Combinatorics
using Optimization
using OptimizationOptimJL
using NonlinearOptimizationTestFunctions
using LinearAlgebra

const NOTF = NonlinearOptimizationTestFunctions
const MIN_EPS = 1e-12

# Exportiere Helper für die Tests
export analyze, SHGOResult, MinimumPoint
export PointCache, KuhnTopology, get_simplices_in_cube, index_to_position, is_valid_index, num_evaluated, get_neighbors, get_value!

# =============================================================================
# Result Types
# =============================================================================

struct MinimumPoint
    minimizer::Vector{Float64}
    objective::Float64
end

Base.getproperty(m::MinimumPoint, s::Symbol) =
    s === :u ? getfield(m, :minimizer) : getfield(m, s)

struct SHGOResult
    local_minima::Vector{MinimumPoint}
    num_basins::Int
    iterations::Int
    converged::Bool
    f_calls::Int
end

Base.getproperty(r::SHGOResult, s::Symbol) =
    s === :results ? getfield(r, :local_minima) : getfield(r, s)

# =============================================================================
# PointCache
# =============================================================================

struct PointCache{N}
    values::Dict{NTuple{N,Int},Float64}
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
        Dict{NTuple{N,Int},Float64}(),
        SVector{N}(lb),
        SVector{N}(ub),
        SVector{N}(divisions),
        step,
        Ref(0)
    )
end

@inline index_to_position(cache::PointCache{N}, idx::NTuple{N,Int}) where N =
    cache.lb .+ SVector{N}(idx) .* cache.step

@inline is_valid_index(cache::PointCache{N}, idx::NTuple{N,Int}) where N =
    @inbounds all(i -> 0 ≤ idx[i] ≤ cache.divisions[i], 1:N)

@inline function get_value!(cache::PointCache{N}, idx::NTuple{N,Int}, f) where N
    # KORREKTUR: Nutzung von is_valid_index statt manueller Schleife
    !is_valid_index(cache, idx) && return Inf
    
    get!(cache.values, idx) do
        cache.eval_count[] += 1
        f(index_to_position(cache, idx))
    end
end

num_evaluated(cache::PointCache) = cache.eval_count[]

# =============================================================================
# Kuhn Topology
# =============================================================================

struct KuhnTopology{N}
    divisions::SVector{N,Int}
end

KuhnTopology(divisions::Vector{Int}) =
    KuhnTopology{length(divisions)}(SVector{length(divisions)}(divisions))

@inline function get_neighbors(::KuhnTopology{N}, idx::NTuple{N,Int}) where N
    neighbors = NTuple{N,Int}[]
    for delta in Iterators.product(ntuple(_ -> (-1,0,1), N)...)
        all(d -> d == 0, delta) && continue
        push!(neighbors, ntuple(i -> idx[i] + delta[i], N))
    end
    neighbors
end

function get_simplices_in_cube(::KuhnTopology{N}, corner::NTuple{N,Int}) where N
    simplices = Vector{Vector{NTuple{N,Int}}}()
    base = collect(corner)
    for p in permutations(1:N)
        simplex = NTuple{N,Int}[]
        current = copy(base)
        push!(simplex, Tuple(current))
        for d in p
            current[d] += 1
            push!(simplex, Tuple(current))
        end
        push!(simplices, simplex)
    end
    simplices
end

# =============================================================================
# Star-Minimum Detection
# =============================================================================

function is_star_minimum(cache::PointCache{N}, topo::KuhnTopology{N},
                         idx::NTuple{N,Int}, f;
                         rel_tol::Float64 = 1e-10) where N
    val = get_value!(cache, idx, f)
    !isfinite(val) && return false
    tol = max(MIN_EPS, abs(val) * rel_tol)
    for nb in get_neighbors(topo, idx)
        get_value!(cache, nb, f) < val - tol && return false
    end
    true
end

function find_star_minima(cache::PointCache{N}, topo::KuhnTopology{N}, f) where N
    minima = NTuple{N,Int}[]
    for idx in Iterators.product((0:d for d in topo.divisions)...)
        t = NTuple{N,Int}(idx)
        is_star_minimum(cache, topo, t, f) && push!(minima, t)
    end
    minima
end

# =============================================================================
# Basin Clustering
# =============================================================================

function cluster_basins(cache::PointCache{N}, topo::KuhnTopology{N},
                        star_minima::Vector{NTuple{N,Int}}, f;
                        threshold_ratio::Float64 = 0.1) where N
    isempty(star_minima) && return Vector{Vector{NTuple{N,Int}}}()
    length(star_minima) == 1 && return [star_minima]

    vals = [v for v in values(cache.values) if isfinite(v)]
    isempty(vals) && return [[m] for m in star_minima]
    value_range = max(maximum(vals) - minimum(vals), MIN_EPS)

    parent = Dict(m => m for m in star_minima)
    rank   = Dict(m => 0 for m in star_minima)

    find_root(x) = parent[x] == x ? x : (parent[x] = find_root(parent[x]))

    function union!(x, y)
        rx, ry = find_root(x), find_root(y)
        rx == ry && return
        if rank[rx] < rank[ry]
            parent[rx] = ry
        elseif rank[rx] > rank[ry]
            parent[ry] = rx
        else
            parent[ry] = rx
            rank[rx] += 1
        end
    end

    star_set = Set(star_minima)
    for m in star_minima
        v = get_value!(cache, m, f)
        for nb in get_neighbors(topo, m)
            nb ∈ star_set || continue
            abs(v - get_value!(cache, nb, f)) < value_range * threshold_ratio &&
                union!(m, nb)
        end
    end

    clusters = Dict{NTuple{N,Int},Vector{NTuple{N,Int}}}()
    for m in star_minima
        push!(get!(clusters, find_root(m), NTuple{N,Int}[]), m)
    end
    collect(values(clusters))
end

# =============================================================================
# Local Optimization
# =============================================================================

function local_optimize(f, grad, x0, lb, ub; maxiters::Int=500)
    eps = 1e-8
    x0s = clamp.(x0, lb .+ eps, ub .- eps)

    optf = OptimizationFunction(
        (x,p) -> f(x);
        grad = (G,x,p) -> copyto!(G, grad(x))
    )
    prob = OptimizationProblem(optf, x0s; lb=lb, ub=ub)

    try
        sol = solve(prob, LBFGS(); maxiters=maxiters)
        return MinimumPoint(Vector(sol.u), sol.objective)
    catch
        try
            sol = solve(prob, BFGS(); maxiters=maxiters)
            return MinimumPoint(Vector(sol.u), sol.objective)
        catch
            sol = solve(prob, NelderMead(); maxiters=maxiters)
            return MinimumPoint(Vector(sol.u), sol.objective)
        end
    end
end

function deduplicate_minima(minima::Vector{MinimumPoint}; 
                           dist_tol::Float64 = 0.05,
                           val_tol::Float64 = 1e-4)
    unique_minima = MinimumPoint[]
    sorted = sort(minima, by = m -> m.objective)
    
    for m in sorted
        is_new = true
        for existing in unique_minima
            if norm(m.minimizer - existing.minimizer) < dist_tol
                is_new = false
                break
            end
        end
        is_new && push!(unique_minima, m)
    end
    return unique_minima
end

# =============================================================================
# Main Entry Point
# =============================================================================

function analyze(tf;
    n_div_initial::Int = 8,
    n_div_max::Int = 25,
    stability_count::Int = 2,
    threshold_ratio::Float64 = 0.1,
    min_distance_tolerance::Float64 = 0.05,
    local_maxiters::Int = 500,
    verbose::Bool = false,
    n_div::Union{Int,Nothing} = nothing,
    use_gradient_pruning::Bool = false,
    kwargs...
)
    !isnothing(n_div) && (n_div_initial = n_div;
                          n_div_max = max(n_div_max, n_div + 10))

    lb = Vector(NOTF.lb(tf))
    ub = Vector(NOTF.ub(tf))
    N = length(lb)

    f    = x -> tf.f(x)
    grad = x -> tf.grad(x)

    prev_basins = -1
    stable = 0
    iteration = 0
    current = n_div_initial

    final_cache = nothing
    final_basins = Vector{Vector{NTuple{N,Int}}}()

    while current ≤ n_div_max
        iteration += 1
        cache = PointCache(lb, ub, fill(current, N))
        topo  = KuhnTopology(fill(current, N))

        stars  = find_star_minima(cache, topo, f)
        basins = cluster_basins(cache, topo, stars, f;
                                threshold_ratio=threshold_ratio)

        stable = length(basins) == prev_basins && !isempty(basins) ? stable + 1 : 0
        prev_basins = length(basins)

        final_cache = cache
        final_basins = basins
        stable ≥ stability_count && break
        current += 2
    end

    candidates = MinimumPoint[]
    for basin in final_basins
        best = basin[argmin(get_value!(final_cache, i, f) for i in basin)]
        x0 = Vector(index_to_position(final_cache, best))
        push!(candidates, local_optimize(f, grad, x0, lb, ub;
                                     maxiters=local_maxiters))
    end

    unique_minima = deduplicate_minima(candidates; dist_tol=min_distance_tolerance)

    SHGOResult(unique_minima, length(unique_minima), iteration,
               stable ≥ stability_count, num_evaluated(final_cache))
end

end # module