module SHGO

using StaticArrays
using Combinatorics
using Optim
using NonlinearOptimizationTestFunctions
using LinearAlgebra
using FunctionWrappers: FunctionWrapper

const NOTF = NonlinearOptimizationTestFunctions
const MIN_EPS = 1e-12

# Fixe Funktionstypen — eliminiert Neukompilierung bei jedem Aufruf
const ObjFunc = FunctionWrapper{Float64, Tuple{Vector{Float64}}}
const GradFunc = FunctionWrapper{Vector{Float64}, Tuple{Vector{Float64}}}

export analyze, SHGOResult, MinimumPoint, ModalityResult
export PointCache, GlobalCache, KuhnTopology, index_to_position, is_valid_index,
       num_evaluated, get_direct_neighbors, get_neighbors, get_value!,
       get_simplices_in_cube

# =============================================================================
# Result Types
# =============================================================================

struct MinimumPoint
    minimizer::Vector{Float64}
    objective::Float64
end

Base.getproperty(m::MinimumPoint, s::Symbol) =
    s === :u ? getfield(m, :minimizer) : getfield(m, s)

struct ModalityResult
    num_basins::Int
    is_unimodal::Bool
    confidence::Symbol  # :high, :medium, :low
    evidence::String
end

struct SHGOResult
    local_minima::Vector{MinimumPoint}
    num_basins::Int
    iterations::Int
    converged::Bool
    f_calls::Int
    global_minimum::Union{MinimumPoint, Nothing}
    modality::Union{ModalityResult, Nothing}
    convergence_history::Vector{Int}
end

# Abwärtskompatibel: 5-Argument-Konstruktor
function SHGOResult(local_minima::Vector{MinimumPoint}, num_basins::Int,
                    iterations::Int, converged::Bool, f_calls::Int)
    gmin = isempty(local_minima) ? nothing : local_minima[1]
    SHGOResult(local_minima, num_basins, iterations, converged, f_calls,
               gmin, nothing, Int[])
end

Base.getproperty(r::SHGOResult, s::Symbol) =
    s === :results ? getfield(r, :local_minima) : getfield(r, s)

# =============================================================================
# PointCache (behalten für Kompatibilität und Tests)
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
    !is_valid_index(cache, idx) && return Inf
    get!(cache.values, idx) do
        cache.eval_count[] += 1
        f(index_to_position(cache, idx))
    end
end

num_evaluated(cache::PointCache) = cache.eval_count[]

# =============================================================================
# GlobalCache mit Rational Keys
# =============================================================================

struct GlobalCache{N}
    values::Dict{NTuple{N,Rational{Int64}}, Float64}
    lb::SVector{N,Float64}
    ub::SVector{N,Float64}
    eval_count::Base.RefValue{Int}
end

function GlobalCache(lb::Vector{Float64}, ub::Vector{Float64})
    N = length(lb)
    GlobalCache{N}(
        Dict{NTuple{N,Rational{Int64}}, Float64}(),
        SVector{N}(lb),
        SVector{N}(ub),
        Ref(0)
    )
end

@inline function normalized_key(idx::NTuple{N,Int}, divisions::SVector{N,Int}) where N
    ntuple(i -> Rational{Int64}(idx[i], divisions[i]), N)
end

@inline function position_from_key(gc::GlobalCache{N}, key::NTuple{N,Rational{Int64}}) where N
    sv = ntuple(i -> Float64(key[i]), N)
    SVector{N,Float64}(sv) .* (gc.ub .- gc.lb) .+ gc.lb
end

@inline function index_to_position(gc::GlobalCache{N}, idx::NTuple{N,Int},
                                   divisions::SVector{N,Int}) where N
    position_from_key(gc, normalized_key(idx, divisions))
end

@inline function get_value!(gc::GlobalCache{N}, key::NTuple{N,Rational{Int64}}, f) where N
    get!(gc.values, key) do
        gc.eval_count[] += 1
        pos = position_from_key(gc, key)
        f(Vector{Float64}(pos))
    end
end

@inline function get_value!(gc::GlobalCache{N}, idx::NTuple{N,Int},
                            divisions::SVector{N,Int}, f) where N
    get_value!(gc, normalized_key(idx, divisions), f)
end

num_evaluated(gc::GlobalCache) = gc.eval_count[]

# =============================================================================
# Kuhn Topology
# =============================================================================

struct KuhnTopology{N}
    divisions::SVector{N,Int}
end

KuhnTopology(divisions::Vector{Int}) =
    KuhnTopology{length(divisions)}(SVector{length(divisions)}(divisions))

@inline function is_valid_index(topo::KuhnTopology{N}, idx::NTuple{N,Int}) where N
    @inbounds all(i -> 0 ≤ idx[i] ≤ topo.divisions[i], 1:N)
end

# Volle Nachbarschaft (3^N - 1)
@inline function get_neighbors(::KuhnTopology{N}, idx::NTuple{N,Int}) where N
    neighbors = NTuple{N,Int}[]
    for delta in Iterators.product(ntuple(_ -> (-1,0,1), N)...)
        all(d -> d == 0, delta) && continue
        push!(neighbors, ntuple(i -> idx[i] + delta[i], N))
    end
    neighbors
end

# Achsenparallele Nachbarn (2N)
@inline function get_direct_neighbors(::KuhnTopology{N}, idx::NTuple{N,Int}) where N
    neighbors = NTuple{N,Int}[]
    sizehint!(neighbors, 2N)
    for d in 1:N
        push!(neighbors, ntuple(i -> i == d ? idx[i] - 1 : idx[i], N))
        push!(neighbors, ntuple(i -> i == d ? idx[i] + 1 : idx[i], N))
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

function is_star_minimum(gc::GlobalCache{N}, topo::KuhnTopology{N},
                         idx::NTuple{N,Int}, f;
                         rel_tol::Float64 = 1e-10) where N
    divs = topo.divisions
    val = get_value!(gc, idx, divs, f)
    !isfinite(val) && return false
    tol = max(MIN_EPS, abs(val) * rel_tol)
    for nb in get_neighbors(topo, idx)
        is_valid_index(topo, nb) || continue
        get_value!(gc, nb, divs, f) < val - tol && return false
    end
    true
end

function is_star_minimum(gc::GlobalCache{N}, topo::KuhnTopology{N},
                         idx::NTuple{N,Int}, f, grad::GradFunc,
                         use_gradient_pruning::Bool;
                         rel_tol::Float64 = 1e-10,
                         grad_tol::Float64 = 0.1) where N
    is_star_minimum(gc, topo, idx, f; rel_tol) || return false
    if use_gradient_pruning
        pos = index_to_position(gc, idx, topo.divisions)
        val = get_value!(gc, idx, topo.divisions, f)
        g = grad(Vector{Float64}(pos))
        norm(g) >= grad_tol * (1 + abs(val)) && return false
    end
    true
end

function find_star_minima(gc::GlobalCache{N}, topo::KuhnTopology{N}, f) where N
    minima = NTuple{N,Int}[]
    for idx in Iterators.product((0:d for d in topo.divisions)...)
        t = NTuple{N,Int}(idx)
        is_star_minimum(gc, topo, t, f) && push!(minima, t)
    end
    minima
end

function find_star_minima(gc::GlobalCache{N}, topo::KuhnTopology{N}, f,
                          grad::GradFunc, use_gradient_pruning::Bool) where N
    !use_gradient_pruning && return find_star_minima(gc, topo, f)
    minima = NTuple{N,Int}[]
    for idx in Iterators.product((0:d for d in topo.divisions)...)
        t = NTuple{N,Int}(idx)
        is_star_minimum(gc, topo, t, f, grad, true) && push!(minima, t)
    end
    minima
end

# =============================================================================
# Basin Clustering (Descent-basiert + Wertbasiert)
# =============================================================================

function assign_basin_by_descent(gc::GlobalCache{N}, topo::KuhnTopology{N},
                                 idx::NTuple{N,Int}, f) where N
    current = idx
    divs = topo.divisions
    max_steps = sum(topo.divisions) * 2

    for _ in 1:max_steps
        val = get_value!(gc, current, divs, f)
        best_nb = current
        best_val = val
        for nb in get_direct_neighbors(topo, current)
            is_valid_index(topo, nb) || continue
            nb_val = get_value!(gc, nb, divs, f)
            if nb_val < best_val
                best_nb = nb
                best_val = nb_val
            end
        end
        best_nb == current && break
        current = best_nb
    end
    return current
end

function cluster_basins(gc::GlobalCache{N}, topo::KuhnTopology{N},
                        star_minima::Vector{NTuple{N,Int}}, f;
                        threshold_ratio::Float64 = 0.1) where N
    isempty(star_minima) && return Vector{Vector{NTuple{N,Int}}}()
    length(star_minima) == 1 && return [star_minima]

    star_set = Set(star_minima)
    divs = topo.divisions

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

    # Schritt 1: Descent von Nicht-Star-Nachbarn
    for m in star_minima
        for nb in get_neighbors(topo, m)
            is_valid_index(topo, nb) || continue
            nb ∈ star_set && continue
            sink = assign_basin_by_descent(gc, topo, nb, f)
            if sink ∈ star_set && sink != m
                union!(m, sink)
            end
        end
    end

    # Schritt 2: Lokales Wert-Mergen (nur direkte Nachbarn)
    all_vals = [get_value!(gc, m, divs, f) for m in star_minima]
    finite_vals = filter(isfinite, all_vals)
    if !isempty(finite_vals)
        value_range = max(maximum(finite_vals) - minimum(finite_vals), MIN_EPS)
        for m in star_minima
            vm = get_value!(gc, m, divs, f)
            for nb in get_neighbors(topo, m)
                nb ∈ star_set || continue
                abs(vm - get_value!(gc, nb, divs, f)) < value_range * threshold_ratio &&
                    union!(m, nb)
            end
        end
    end

    clusters = Dict{NTuple{N,Int}, Vector{NTuple{N,Int}}}()
    for m in star_minima
        push!(get!(clusters, find_root(m), NTuple{N,Int}[]), m)
    end
    collect(values(clusters))
end

# =============================================================================
# Local Optimization — Zwei-Phasen-Architektur
#
#   Phase "explore": maxiters=100, g_tol=1e-6  → schnelle Basinerkundung
#   Phase "refine":  maxiters=1000, g_tol=1e-10 → präzise Top-Kandidaten
#
# Warum: mishra3 hat 28 Basins. Vorher 28×500 = 14000 L-BFGS-Iterationen.
# Jetzt: explore(28×100) + refine(5×1000) = 2800 + 5000 = 7800.
# Gleiche Genauigkeit, halbe Zeit.
# =============================================================================

function _optimize_lbfgs(f::F, g!::G, lb, ub, x0, maxiters;
                         g_tol::Float64 = 1e-6) where {F,G}
    optimize(f, g!, lb, ub, x0, Fminbox(LBFGS()),
             Optim.Options(iterations=maxiters, g_tol=g_tol,
                           allow_f_increases=true))
end

function _optimize_nm(f::F, lb, ub, x0, maxiters) where {F}
    optimize(f, lb, ub, x0, Fminbox(NelderMead()),
             Optim.Options(iterations=min(maxiters, 500)))
end

const _dummy_f = x -> sum(x.^2)
const _dummy_g! = (G, x) -> (G .= 2 .* x)
const _precompile_done = Ref(false)

function _ensure_precompiled()
    _precompile_done[] && return
    try
        _optimize_lbfgs(_dummy_f, _dummy_g!, [0.0, 0.0], [1.0, 1.0], [0.5, 0.5], 10)
    catch; end
    _precompile_done[] = true
end

"""
    local_optimize(f, grad, x0, lb, ub; maxiters=500, g_tol=1e-8)

L-BFGS mit Fminbox. NelderMead als Fallback bei Fehler.
"""
function local_optimize(f, grad, x0::Vector{Float64}, lb::Vector{Float64},
                        ub::Vector{Float64};
                        maxiters::Int = 500, g_tol::Float64 = 1e-8)
    _ensure_precompiled()
    x0_clamped = clamp.(x0, lb .+ 1e-10, ub .- 1e-10)
    g! = (G, x) -> (G .= grad(x))

    try
        result = _optimize_lbfgs(f, g!, lb, ub, x0_clamped, maxiters; g_tol)
        return MinimumPoint(Optim.minimizer(result), Optim.minimum(result))
    catch
        try
            result = _optimize_nm(f, lb, ub, x0_clamped, maxiters)
            return MinimumPoint(Optim.minimizer(result), Optim.minimum(result))
        catch
            return MinimumPoint(x0_clamped, f(x0_clamped))
        end
    end
end

"""
    optimize_basin(basin, gc, topo, f, grad, lb, ub; k, maxiters, g_tol)

Multi-Start pro Basin. Gibt ALLE k Minima zurück (nicht nur das beste),
damit deduplicate_minima korrekt zählen kann.
"""
function optimize_basin(basin::Vector{NTuple{N,Int}}, gc::GlobalCache{N},
                        topo::KuhnTopology{N}, f, grad,
                        lb::Vector{Float64}, ub::Vector{Float64};
                        k::Int = 3, maxiters::Int = 500,
                        g_tol::Float64 = 1e-8) where N
    divs = topo.divisions
    sorted = sort(basin, by = idx -> get_value!(gc, idx, divs, f))
    n_starts = min(k, length(sorted))

    results = MinimumPoint[]
    for i in 1:n_starts
        x0 = Vector{Float64}(index_to_position(gc, sorted[i], divs))
        push!(results, local_optimize(f, grad, x0, lb, ub; maxiters, g_tol))
    end
    results
end

"""
    generate_probe_points(lb, ub) → Vector{Vector{Float64}}

Strategische Startpunkte die das Kuhn-Grid ergänzen:
  - Origin (0,...,0) falls innerhalb Bounds → bohachevsky, step, schwefel
  - Center of bounds → symmetrische Funktionen
  - Quartil-Punkte entlang jeder Achse → abseits Grid-Linien
  - Extreme Ecken → Randbereiche des Suchraums

Kosten: 2N+4 Punkte (N=2: 8, N=3: 10, N=4: 12)
"""
function generate_probe_points(lb::Vector{Float64}, ub::Vector{Float64})
    N = length(lb)
    probes = Vector{Float64}[]

    # Origin (0,...,0) — wenn innerhalb Bounds
    origin = zeros(N)
    if all(i -> lb[i] ≤ 0.0 ≤ ub[i], 1:N)
        push!(probes, origin)
    end

    # Center
    center = (lb .+ ub) ./ 2
    push!(probes, center)

    # Quartil-Punkte: 25% und 75% entlang jeder Achse, Rest auf Center
    for d in 1:N
        for frac in (0.25, 0.75)
            p = copy(center)
            p[d] = lb[d] + frac * (ub[d] - lb[d])
            push!(probes, p)
        end
    end

    # Extreme Ecken (leicht nach innen versetzt)
    eps_b = (ub .- lb) .* 1e-4
    push!(probes, lb .+ eps_b)
    push!(probes, ub .- eps_b)

    probes
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
# Modality Assessment
# =============================================================================

function assess_modality(num_basins::Int, converged::Bool, iterations::Int)
    if num_basins == 1 && converged
        ModalityResult(1, true, :high,
            "Einziges Basin, Betti-Zahl stabil über $iterations Iterationen")
    elseif num_basins == 1 && !converged
        ModalityResult(1, true, :medium,
            "Einziges Basin gefunden, aber Konvergenz nicht bestätigt")
    elseif num_basins > 1 && converged
        ModalityResult(num_basins, false, :high,
            "$num_basins stabile Basins — Funktion ist multimodal")
    else
        ModalityResult(num_basins, false, :low,
            "$num_basins Basins, aber noch nicht stabilisiert")
    end
end

# =============================================================================
# Main Entry Point
#
# Architektur:
#   Phase 1: Topologische Analyse (Grid → Stars → Descent-Clustering → Basins)
#   Phase 2: Schnelle Basinerkundung (maxiters=100, budget-limitiert)
#   Phase 3: Strategische Probes (Origin, Center, Quartile, Ecken)
#   Phase 4: Präzise Verfeinerung (Top-N Kandidaten, maxiters=1000)
#   Phase 5: Deduplizierung → Modalitäts-Assessment → Ergebnis
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
    multi_start_k::Int = 3,
    max_local_opts::Int = 30,
    max_stars::Int = 50,
    n_refine::Int = 5,
    kwargs...
)
    _ensure_precompiled()

    !isnothing(n_div) && (n_div_initial = n_div;
                          n_div_max = max(n_div_max, n_div + 10))

    lb = Float64.(NOTF.lb(tf))
    ub = Float64.(NOTF.ub(tf))
    N = length(lb)

    f = ObjFunc(x -> Float64(tf.f(x)))
    grad = GradFunc(x -> Vector{Float64}(tf.grad(x)))

    gc = GlobalCache(lb, ub)

    # =====================================================================
    # PHASE 1: Topologische Analyse
    # =====================================================================
    prev_basins = -1
    stable = 0
    iteration = 0
    current = n_div_initial
    convergence_history = Int[]

    final_basins = Vector{Vector{NTuple{N,Int}}}()
    final_topo = nothing

    while current ≤ n_div_max
        iteration += 1
        topo = KuhnTopology(fill(current, N))

        stars = find_star_minima(gc, topo, f, grad, use_gradient_pruning)

        # Star-Limit: verhindert Clustering-Explosion
        if length(stars) > max_stars
            divs = topo.divisions
            sort!(stars, by = idx -> get_value!(gc, idx, divs, f))
            resize!(stars, max_stars)
        end

        basins = cluster_basins(gc, topo, stars, f;
                                threshold_ratio=threshold_ratio)

        num_b = length(basins)
        push!(convergence_history, num_b)

        if verbose
            println("  Iteration $iteration: n_div=$current, " *
                    "stars=$(length(stars)), basins=$num_b, " *
                    "evals=$(num_evaluated(gc))")
        end

        stable = num_b == prev_basins && num_b > 0 ? stable + 1 : 0
        prev_basins = num_b

        final_basins = basins
        final_topo = topo
        stable ≥ stability_count && break

        # Early exit: hochoszillierende Funktionen
        if iteration ≥ 3 && num_b > 15 && stable == 0
            verbose && println("  Early exit: $num_b unstable basins after $iteration iters")
            break
        end

        current += 2
    end

    converged = stable ≥ stability_count

    # =====================================================================
    # PHASE 2: Schnelle Basinerkundung (maxiters=100)
    # =====================================================================
    candidates = MinimumPoint[]

    if final_topo !== nothing && !isempty(final_basins)
        divs = final_topo.divisions
        n_basins = length(final_basins)

        basin_scores = [(i, minimum(get_value!(gc, idx, divs, f) for idx in basin))
                        for (i, basin) in enumerate(final_basins)]
        sort!(basin_scores, by = x -> x[2])

        budget = max_local_opts
        k_base = max(1, min(multi_start_k, budget ÷ max(n_basins, 1)))

        for (i, _) in basin_scores
            budget <= 0 && break
            k = min(k_base, budget, length(final_basins[i]))
            k = max(1, k)
            basin_results = optimize_basin(final_basins[i], gc, final_topo,
                                           f, grad, lb, ub;
                                           k=k,
                                           maxiters=100,
                                           g_tol=1e-6)
            append!(candidates, basin_results)
            budget -= k
        end

        verbose && println("  Phase 2: $(max_local_opts - budget)/$max_local_opts starts, $n_basins basins")
    end

    # =====================================================================
    # PHASE 3: Strategische Probes
    # =====================================================================
    probes = generate_probe_points(lb, ub)
    for x0 in probes
        push!(candidates, local_optimize(f, grad, x0, lb, ub;
                                         maxiters=200, g_tol=1e-8))
    end

    # Zusätzliche Random-Starts bei fehlender Konvergenz
    if !converged
        n_extra = max(5, 3 * N)
        for _ in 1:n_extra
            x0 = lb .+ rand(N) .* (ub .- lb)
            push!(candidates, local_optimize(f, grad, x0, lb, ub;
                                             maxiters=200, g_tol=1e-8))
        end
    end

    # =====================================================================
    # PHASE 4: Präzise Verfeinerung der Top-Kandidaten
    # =====================================================================
    diag = norm(ub .- lb)
    dyn_dist_tol = max(min_distance_tolerance, diag * 0.001)
    unique_minima = deduplicate_minima(candidates; dist_tol=dyn_dist_tol)

    n_ref = min(n_refine, length(unique_minima))
    for i in 1:n_ref
        m = unique_minima[i]
        refined = local_optimize(f, grad, m.minimizer, lb, ub;
                                 maxiters=1000, g_tol=1e-10)
        if refined.objective ≤ m.objective
            unique_minima[i] = refined
        end
    end

    # Re-sort und re-dedup nach Verfeinerung
    sort!(unique_minima, by = m -> m.objective)
    unique_minima = deduplicate_minima(unique_minima; dist_tol=dyn_dist_tol)

    # =====================================================================
    # PHASE 5: Ergebnis
    # =====================================================================
    num_basins_final = length(unique_minima)
    global_min = isempty(unique_minima) ? nothing : unique_minima[1]
    modality = assess_modality(num_basins_final, converged, iteration)

    SHGOResult(unique_minima, num_basins_final, iteration, converged,
               num_evaluated(gc), global_min, modality, convergence_history)
end

function __init__()
    _ensure_precompiled()
end

end # module
