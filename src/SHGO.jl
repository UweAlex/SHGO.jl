module SHGO

using StaticArrays, Combinatorics, Optimization, OptimizationOptimJL, NonlinearOptimizationTestFunctions, LazySets, LinearAlgebra, Printf

# Inkludieren der Teilmodule (Basins wird nicht mehr benötigt)
include("triangulation/Grid.jl")
include("triangulation/TopicalStructure.jl")

const G    = Grid
const TS   = TopicalStructure
const NOTF = NonlinearOptimizationTestFunctions

export analyze, SHGOResult, MinimumPoint

"""
    MinimumPoint
Struktur für einen lokalen Minimierer. Das Feld .objective wird vom Skript benötigt.
"""
struct MinimumPoint
    minimizer::Vector{Float64}
    objective::Float64
end

# Ermöglicht den Zugriff auf .u (SciML Standard) und .minimizer (Skript Standard)
Base.getproperty(m::MinimumPoint, s::Symbol) = s === :u ? getfield(m, :minimizer) : getfield(m, s)

"""
    SHGOResult
Ergebniscontainer. num_basins wird hier über die Anzahl unique Minima definiert.
"""
struct SHGOResult
    results::Vector{MinimumPoint}
    num_basins::Int
    local_minima::Vector{MinimumPoint}
end

mutable struct SHGOState{N, TF, GTYPE <: G.AbstractGrid}
    tf::TF
    grid::GTYPE
    tm::TS.TopicalManager{N}
end

"""
    analyze(tf; n_div, verbose, use_gradient_pruning, refinement_levels, min_distance_tolerance)

Hauptfunktion, korrigiert für Diagnose-Skript-Kompatibilität.
"""
function analyze(tf; 
                 n_div::Int=10, 
                 verbose::Bool=false, 
                 use_gradient_pruning::Bool=true,
                 refinement_levels::Int=0,  # Korrektur: Argument wird jetzt akzeptiert!
                 min_distance_tolerance::Float64=0.01)
    
    if verbose
        @printf "--- SHGO: Starte Analyse (n_div=%d, pruning=%s) ---\n" n_div use_gradient_pruning
    end

    # 1. Setup
    lb = Vector{Float64}(NOTF.lb(tf)); ub = Vector{Float64}(NOTF.ub(tf)); N = length(lb)
    grid = G.GridStructure(lb, ub, fill(n_div, N))
    tm = TS.TopicalManager{N}()
    state = SHGOState{N, typeof(tf), typeof(grid)}(tf, grid, tm)

    # 2. Triangulation
    perms = collect(permutations(1:N))
    ranges = ntuple(i -> 1:(state.grid.dims[i]-1), N)
    
    for cube_idx in CartesianIndices(ranges)
        for p in perms
            v_ids = Int[]; curr = cube_idx
            for step in 0:N
                if haskey(state.grid.cache.points, curr)
                    push!(v_ids, state.grid.cache.points[curr].v_id)
                else
                    pos = G.calculate_pos(state.grid, curr)
                    val = tf.f(pos)
                    v = TS.add_vertex!(state.tm, pos, tf.f(pos))
                    state.grid.cache.points[curr] = G.GridPoint(curr, pos, v.id)
                    push!(v_ids, v.id)
                end
                if step < N
                    delta = ntuple(i -> i == p[step+1] ? 1 : 0, N)
                    curr += CartesianIndex(delta)
                end
            end
            TS.add_simplex!(state.tm, v_ids)
        end
    end

    # 3. Pruning
    active_s_ids = use_gradient_pruning ? prune!(state) : TS.all_simplex_ids(state.tm)
    
    # 4. Minimizer Pool (SciPy-Style)
    candidate_minima = MinimumPoint[]
    for s_id in active_s_ids
        s = state.tm.simplices[s_id]
        # Startpunkt ist der beste Vertex im Simplex
        best_v_id = s.vertices[argmin([state.tm.vertices[v].val for v in s.vertices])]
        x0 = Vector{Float64}(state.tm.vertices[best_v_id].pos)
        
        try
            f_opt = OptimizationFunction((x,p) -> tf.f(x), Optimization.AutoForwardDiff())
            prob = OptimizationProblem(f_opt, x0; lb=lb, ub=ub)
            sol = solve(prob, LBFGS(); maxiters=100)
            push!(candidate_minima, MinimumPoint(Vector{Float64}(sol.minimizer), Float64(sol.objective)))
        catch; end
    end

    # 5. Deduplizierung
    unique_minima = remove_duplicates(candidate_minima, min_distance_tolerance)
    sort!(unique_minima, by = m -> m.objective)
    
    return SHGOResult(unique_minima, length(unique_minima), unique_minima)
end

function remove_duplicates(minima::Vector{MinimumPoint}, tol::Float64)
    isempty(minima) && return MinimumPoint[]
    sorted = sort(minima, by = m -> m.objective)
    unique_list = MinimumPoint[sorted[1]]
    for cand in sorted[2:end]
        if all(norm(cand.minimizer - u.minimizer) > tol for u in unique_list)
            push!(unique_list, cand)
        end
    end
    return unique_list
end

function prune!(state::SHGOState{N}) where N
    active = Int[]
    for s_id in TS.all_simplex_ids(state.tm)
        s = state.tm.simplices[s_id]
        grads = [Vector{Float64}(state.tf.grad(state.tm.vertices[vid].pos)) for vid in s.vertices]
        if zeros(N) ∈ VPolytope(grads)
            push!(active, s_id)
        end
    end
    return active
end

end # module