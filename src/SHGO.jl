# File: src/SHGO.jl
#
# SHGO.jl - Simplicial Homology Global Optimization
#
# Hochperformante Julia-Implementierung basierend auf:
# - Implizite Kuhn-Topologie (kein expliziter Graph im Speicher)
# - Lazy Evaluation mit Point-Cache
# - Betti-Zahl-Stabilität als Konvergenzkriterium
# - Star-Domain basierte Minimum-Kandidaten-Erkennung
# - Infinity-Padding für Randbehandlung
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
# Ergebnis-Typen
# =============================================================================

struct MinimumPoint
    minimizer::Vector{Float64}
    objective::Float64
end

Base.getproperty(m::MinimumPoint, s::Symbol) =
    s === :u ? getfield(m, :minimizer) : getfield(m, s)

struct SHGOResult
    results::Vector{MinimumPoint}
    num_basins::Int
    local_minima::Vector{MinimumPoint}
    iterations::Int
    converged::Bool
end

# =============================================================================
# Point Cache - Effiziente Speicherung evaluierter Punkte
# =============================================================================

"""
    PointCache{N}

Speichert evaluierte Funktionswerte indiziert durch Gitter-Koordinaten.
Verhindert redundante Evaluationen bei überlappenden Simplizes.
"""
struct PointCache{N}
    values::Dict{NTuple{N,Int}, Float64}
    lb::SVector{N,Float64}
    ub::SVector{N,Float64}
    divisions::SVector{N,Int}
    step::SVector{N,Float64}
end

function PointCache(lb::Vector{Float64}, ub::Vector{Float64}, divisions::Vector{Int})
    N = length(lb)
    step = SVector{N}((ub .- lb) ./ divisions)
    PointCache{N}(
        Dict{NTuple{N,Int}, Float64}(),
        SVector{N}(lb),
        SVector{N}(ub),
        SVector{N}(divisions),
        step
    )
end

"""
Konvertiert Gitter-Index zu physikalischer Position.
"""
function index_to_position(cache::PointCache{N}, idx::NTuple{N,Int}) where N
    cache.lb .+ SVector{N}(idx) .* cache.step
end

"""
Prüft ob Index innerhalb der Bounds liegt.
"""
function is_valid_index(cache::PointCache{N}, idx::NTuple{N,Int}) where N
    all(i -> 0 <= idx[i] <= cache.divisions[i], 1:N)
end

"""
Holt oder berechnet Funktionswert. Gibt Inf für Out-of-Bounds (Infinity-Padding).
"""
function get_value!(cache::PointCache{N}, idx::NTuple{N,Int}, f::Function) where N
    # Infinity-Padding: Außerhalb der Bounds → +Inf
    if !is_valid_index(cache, idx)
        return Inf
    end
    
    # Cache-Lookup
    get!(cache.values, idx) do
        pos = index_to_position(cache, idx)
        f(pos)
    end
end

"""
Anzahl der bereits evaluierten Punkte.
"""
num_evaluated(cache::PointCache) = length(cache.values)

# =============================================================================
# Implizite Kuhn-Topologie
# =============================================================================

"""
    KuhnTopology{N}

Implizite Repräsentation der Kuhn-Triangulation.
Berechnet Nachbarschaften on-the-fly ohne expliziten Graphen.
"""
struct KuhnTopology{N}
    divisions::SVector{N,Int}
    # Vorberechnete Permutationen für Kuhn-Simplizes
    permutations::Vector{Vector{Int}}
end

function KuhnTopology(divisions::Vector{Int})
    N = length(divisions)
    perms = collect(permutations(1:N))
    KuhnTopology{N}(SVector{N}(divisions), perms)
end

"""
Gibt alle Kuhn-Nachbarn eines Gitterpunkts zurück.
Ein Nachbar ist ein Punkt, der mit diesem Punkt in mindestens einem Kuhn-Simplex liegt.

Die Kuhn-Triangulation verbindet Punkt (i₁,...,iₙ) mit:
- Allen Punkten die sich in genau einer Koordinate um ±1 unterscheiden
- Plus diagonale Verbindungen entlang der Kuhn-Pfade
"""
function get_neighbors(topo::KuhnTopology{N}, idx::NTuple{N,Int}) where N
    neighbors = Set{NTuple{N,Int}}()
    
    # Achsen-Nachbarn (±1 in jeder Dimension)
    for d in 1:N
        for delta in (-1, 1)
            neighbor = ntuple(i -> i == d ? idx[i] + delta : idx[i], N)
            push!(neighbors, neighbor)
        end
    end
    
    # Diagonale Nachbarn (Kuhn-Pfad-Verbindungen)
    # In der Kuhn-Triangulation sind auch bestimmte Diagonalen verbunden
    for d1 in 1:N, d2 in (d1+1):N
        # (+1, +1) Diagonale
        push!(neighbors, ntuple(i -> i == d1 || i == d2 ? idx[i] + 1 : idx[i], N))
        # (-1, -1) Diagonale
        push!(neighbors, ntuple(i -> i == d1 || i == d2 ? idx[i] - 1 : idx[i], N))
    end
    
    return neighbors
end

"""
Generiert alle Kuhn-Simplizes in einem Hyperwürfel mit Ecke bei `cube_idx`.
Jeder Simplex wird als Tupel von N+1 Vertex-Indizes zurückgegeben.
"""
function get_simplices_in_cube(topo::KuhnTopology{N}, cube_idx::NTuple{N,Int}) where N
    simplices = Vector{NTuple{N+1, NTuple{N,Int}}}()
    
    for perm in topo.permutations
        vertices = Vector{NTuple{N,Int}}(undef, N+1)
        current = cube_idx
        vertices[1] = current
        
        for (step, dim) in enumerate(perm)
            current = ntuple(i -> i == dim ? current[i] + 1 : current[i], N)
            vertices[step + 1] = current
        end
        
        push!(simplices, NTuple{N+1, NTuple{N,Int}}(Tuple(vertices)))
    end
    
    return simplices
end

# =============================================================================
# Star-Domain basierte Minimum-Erkennung
# =============================================================================

"""
    is_star_minimum(cache, topo, idx, f)

Prüft ob ein Punkt das Minimum in seiner Star-Domain ist.
Die Star-Domain sind alle Simplizes, die diesen Punkt enthalten.

Ein Punkt ist ein lokaler Minimum-Kandidat, wenn sein Funktionswert
kleiner oder gleich allen seinen Kuhn-Nachbarn ist.
"""
function is_star_minimum(
    cache::PointCache{N}, 
    topo::KuhnTopology{N}, 
    idx::NTuple{N,Int},
    f::Function
) where N
    # Eigener Wert
    val = get_value!(cache, idx, f)
    
    # Vergleiche mit allen Nachbarn
    for neighbor in get_neighbors(topo, idx)
        neighbor_val = get_value!(cache, neighbor, f)
        if neighbor_val < val - MIN_EPS
            return false
        end
    end
    
    return true
end

"""
Findet alle Star-Minima im aktuellen Gitter.
"""
function find_star_minima(
    cache::PointCache{N},
    topo::KuhnTopology{N},
    f::Function
) where N
    minima = NTuple{N,Int}[]
    
    # Iteriere über alle Gitterpunkte
    for idx in Iterators.product((0:d for d in topo.divisions)...)
        idx_tuple = NTuple{N,Int}(idx)
        if is_star_minimum(cache, topo, idx_tuple, f)
            push!(minima, idx_tuple)
        end
    end
    
    return minima
end

# =============================================================================
# Basin-Clustering - Direkte Nachbarschaft statt BFS
# =============================================================================

"""
    cluster_basins(cache, topo, star_minima, f; merge_tolerance=0.01)

Clustert Star-Minima zu Basins.

NEUE STRATEGIE:
- Jedes Star-Minimum startet als eigenes Basin
- Zwei Star-Minima werden NUR verbunden, wenn:
  1. Sie direkte Kuhn-Nachbarn sind, ODER
  2. Sie durch eine "flache" Region verbunden sind (Barriere < merge_tolerance * range)

Das verhindert, dass separate Täler fälschlich zusammengeführt werden.
"""
function cluster_basins(
    cache::PointCache{N},
    topo::KuhnTopology{N},
    star_minima::Vector{NTuple{N,Int}},
    f::Function;
    threshold_ratio::Float64 = 0.1  # Legacy-Parameter, wird für merge_tolerance genutzt
) where N
    isempty(star_minima) && return Vector{Vector{NTuple{N,Int}}}()
    
    # Nur 1 Star-Minimum → 1 Basin
    if length(star_minima) == 1
        return [star_minima]
    end
    
    # Berechne Merge-Toleranz basierend auf lokalen Unterschieden
    all_vals = collect(values(cache.values))
    f_min = minimum(all_vals)
    f_max = maximum(all_vals)
    value_range = f_max - f_min
    
    # Union-Find
    parent = Dict{NTuple{N,Int}, NTuple{N,Int}}()
    for m in star_minima
        parent[m] = m
    end
    
    function find_root(x)
        if parent[x] != x
            parent[x] = find_root(parent[x])
        end
        return parent[x]
    end
    
    function union!(x, y)
        rx, ry = find_root(x), find_root(y)
        if rx != ry
            parent[rx] = ry
        end
    end
    
    star_set = Set(star_minima)
    
    # Strategie: Verbinde nur DIREKTE Nachbarn unter bestimmten Bedingungen
    for m1 in star_minima
        val1 = get_value!(cache, m1, f)
        neighbors_of_m1 = get_neighbors(topo, m1)
        
        for m2 in star_minima
            m1 >= m2 && continue  # Nur einmal pro Paar
            
            val2 = get_value!(cache, m2, f)
            
            # Fall 1: Direkte Kuhn-Nachbarn
            if m2 in neighbors_of_m1
                # Verbinde nur wenn die Werte sehr ähnlich sind
                # (= sie gehören zum selben flachen Tal)
                val_diff = abs(val1 - val2)
                if val_diff < value_range * threshold_ratio
                    union!(m1, m2)
                end
                continue
            end
            
            # Fall 2: Prüfe ob ein gemeinsamer Nachbar existiert, der niedriger ist als beide
            # Das würde bedeuten, sie teilen sich ein Tal
            common_neighbors = intersect(neighbors_of_m1, get_neighbors(topo, m2))
            
            for cn in common_neighbors
                cn_val = get_value!(cache, cn, f)
                # Gemeinsamer Nachbar ist niedriger als beide Star-Minima?
                # Das sollte nicht passieren (sonst wären m1/m2 keine Star-Minima)
                # Also: Verbinde nur wenn der gemeinsame Nachbar ähnlich niedrig ist
                if cn_val <= max(val1, val2) + value_range * threshold_ratio * 0.5
                    union!(m1, m2)
                    break
                end
            end
        end
    end
    
    # Gruppiere nach Root
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
# Lokale Optimierung
# =============================================================================

"""
Optimiert lokal von einem Startpunkt aus.
"""
function local_optimize(
    f::Function,
    grad::Function,
    x0::Vector{Float64},
    lb::Vector{Float64},
    ub::Vector{Float64};
    maxiters::Int = 500
)
    fopt = OptimizationFunction(
        (x, p) -> f(x),
        grad = (G, x, p) -> copyto!(G, grad(x))
    )
    
    prob = OptimizationProblem(fopt, x0; lb=lb, ub=ub)
    sol = solve(prob, LBFGS(); maxiters=maxiters)
    
    return MinimumPoint(Vector(sol.minimizer), sol.objective)
end

"""
Dedupliziert Minima basierend auf räumlicher Nähe.
"""
function deduplicate_minima(minima::Vector{MinimumPoint}; dist_tol::Float64 = 0.05)
    isempty(minima) && return MinimumPoint[]
    
    sorted = sort(minima, by = m -> m.objective)
    unique_minima = [sorted[1]]
    
    for m in sorted[2:end]
        is_new = all(u -> norm(m.minimizer - u.minimizer) >= dist_tol, unique_minima)
        is_new && push!(unique_minima, m)
    end
    
    return unique_minima
end

# =============================================================================
# Hauptfunktion: analyze() mit Betti-Stabilität
# =============================================================================

"""
    analyze(tf; kwargs...)

Analysiert die Optimierungslandschaft einer Testfunktion.

# Algorithmus:
1. Starte mit grobem Gitter
2. Finde Star-Minima (lokale Minimum-Kandidaten)
3. Clustere zu Basins
4. Verfeinere Gitter und wiederhole
5. Stoppe wenn Betti-Zahl (Anzahl Basins) stabil ist
6. Optimiere lokal pro Basin

# Keyword Arguments
- `n_div_initial::Int = 8`: Initiale Grid-Auflösung
- `n_div_max::Int = 25`: Maximale Grid-Auflösung
- `stability_count::Int = 2`: Anzahl stabiler Iterationen für Konvergenz
- `threshold_ratio::Float64 = 0.1`: Schwellenwert für Basin-Clustering
- `min_distance_tolerance::Float64 = 0.05`: Deduplikations-Toleranz
- `local_maxiters::Int = 500`: Max. Iterationen für lokale Optimierung
- `verbose::Bool = false`: Debug-Ausgaben

# Legacy Parameters (für Rückwärtskompatibilität)
- `n_div`: Alias für n_div_initial
- `use_gradient_pruning`, `pruning_tol`, etc.: Werden ignoriert
"""
function analyze(
    tf;
    # Neue Parameter
    n_div_initial::Int = 8,
    n_div_max::Int = 25,
    stability_count::Int = 2,
    threshold_ratio::Float64 = 0.1,
    min_distance_tolerance::Float64 = 0.05,
    local_maxiters::Int = 500,
    verbose::Bool = false,
    # Legacy-Parameter (Rückwärtskompatibilität)
    n_div::Union{Int,Nothing} = nothing,
    barrier_tolerance_ratio::Float64 = 0.1,
    use_gradient_pruning::Bool = false,
    pruning_tol::Float64 = 0.2,
    adaptive_refinement::Bool = true,
    adaptive_max_levels::Int = 5,
    refinement_focus::Int = 15,
    triangulation = nothing
)
    # Legacy: n_div überschreibt n_div_initial
    if !isnothing(n_div)
        n_div_initial = n_div
        n_div_max = max(n_div_max, n_div + 10)
    end
    
    # Legacy: barrier_tolerance_ratio → threshold_ratio
    threshold_ratio = barrier_tolerance_ratio
    
    lb = Vector{Float64}(NOTF.lb(tf))
    ub = Vector{Float64}(NOTF.ub(tf))
    N = length(lb)
    
    # Wrapper-Funktionen
    f = x -> tf.f(x)
    grad = x -> tf.grad(x)
    
    # Warnungen für hohe Dimensionen
    if N > 6 && verbose
        @warn "Dimension N=$N ist hoch. Kuhn-Triangulation skaliert mit N!. " *
              "Erwäge Sobol-Sampling für N > 6."
    end
    
    # =========================================================================
    # Iterative Verfeinerung mit Betti-Stabilität
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
        
        # Erstelle Cache und Topologie für aktuelle Auflösung
        divisions = fill(current_n_div, N)
        cache = PointCache(lb, ub, divisions)
        topo = KuhnTopology(divisions)
        
        # Finde Star-Minima
        star_minima = find_star_minima(cache, topo, f)
        
        if verbose
            println("  Star-Minima gefunden: $(length(star_minima))")
            println("  Punkte evaluiert: $(num_evaluated(cache))")
        end
        
        # Clustere zu Basins
        basins = cluster_basins(cache, topo, star_minima, f; 
                               threshold_ratio=threshold_ratio)
        
        num_basins = length(basins)
        
        if verbose
            println("  Basins: $num_basins")
        end
        
        # Prüfe Betti-Stabilität
        if num_basins == prev_num_basins
            stable_iterations += 1
            if verbose
                println("  Stabil seit $stable_iterations Iterationen")
            end
        else
            stable_iterations = 0
        end
        
        prev_num_basins = num_basins
        final_basins = basins
        final_cache = cache
        
        # Konvergenz erreicht?
        if stable_iterations >= stability_count
            if verbose
                println("Konvergenz nach $iteration Iterationen bei n_div=$current_n_div")
            end
            break
        end
        
        # Verfeinere Gitter
        current_n_div += 2
    end
    
    converged = stable_iterations >= stability_count
    
    if !converged && verbose
        @warn "Keine Konvergenz erreicht. Maximale Auflösung n_div=$n_div_max verwendet."
    end
    
    # =========================================================================
    # Lokale Optimierung pro Basin
    # =========================================================================
    
    if verbose
        println("\nLokale Optimierung für $(length(final_basins)) Basins...")
    end
    
    candidates = MinimumPoint[]
    
    for (basin_id, basin) in enumerate(final_basins)
        # Wähle besten Punkt im Basin als Startpunkt
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
                @warn "Lokale Optimierung fehlgeschlagen für Basin $basin_id" exception=e
            end
        end
    end
    
    # =========================================================================
    # Deduplikation
    # =========================================================================
    
    unique_minima = deduplicate_minima(candidates; dist_tol=min_distance_tolerance)
    
    if verbose
        println("\nErgebnis: $(length(unique_minima)) eindeutige Minima")
        for (i, m) in enumerate(sort(unique_minima, by=x->x.objective))
            println("  $i. f = $(round(m.objective, digits=6)) @ $(round.(m.minimizer, digits=4))")
        end
    end
    
    return SHGOResult(
        unique_minima,
        length(unique_minima),
        unique_minima,
        iteration,
        converged
    )
end

end # module