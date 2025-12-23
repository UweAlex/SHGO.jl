# src/triangulation/kuhn.jl
using StaticArrays
using LazySets

# ────────────────────────────────────────────────────────────────
# Heap's Algorithm – Permutations-Generator
# ────────────────────────────────────────────────────────────────

struct KuhnPermutationIterator{N} end

function KuhnPermutationIterator(n::Int)
    n < 0 && throw(ArgumentError("n muss >= 0 sein"))
    KuhnPermutationIterator{n}()
end

# Initialer Aufruf (ohne state)
function Base.iterate(it::KuhnPermutationIterator{N}) where N
    N == 0 && return nothing

    p = MVector{N, Int}(1:N)
    c = MVector{N, Int}(zeros(Int, N))
    return SVector{N, Int}(p), (p, c, 1)
end

# Fortsetzung mit state (nur für gültigen Tuple-Zustand)
function Base.iterate(it::KuhnPermutationIterator{N}, state::Tuple) where N
    p, c, k = state

    while k ≤ N
        if c[k] < k - 1
            if isodd(k)
                p[1], p[k] = p[k], p[1]
            else
                j = c[k] + 1
                p[j], p[k] = p[k], p[j]
            end

            c[k] += 1
            return SVector{N, Int}(p), (p, c, 1)  # Reset k
        else
            c[k] = 0
            k += 1
        end
    end

    nothing
end

Base.length(::KuhnPermutationIterator{N}) where N = factorial(N)
Base.eltype(::Type{KuhnPermutationIterator{N}}) where N = SVector{N, Int}

# ────────────────────────────────────────────────────────────────
# 2. Kuhn-Indizes aus Ursprung + Permutation
# ────────────────────────────────────────────────────────────────

function generate_kuhn_indices(origin::CartesianIndex{N}, perm::SVector{N,Int}) where N
    indices = Vector{CartesianIndex{N}}(undef, N+1)
    indices[1] = origin

    current = origin
    for j in 1:N
        dim = perm[j]
        offset = ntuple(i -> i == dim ? 1 : 0, Val(N))
        current += CartesianIndex(offset)
        indices[j+1] = current
    end

    indices
end

# ────────────────────────────────────────────────────────────────
# 3. Lazy Iterator mit Gradient-Hull-Pruning (deterministisch, stabil)
# ────────────────────────────────────────────────────────────────

struct LazyKuhnSimplexes{N}
    cell_origin_idx::CartesianIndex{N}
    cache::VertexCache{N}
end

Base.IteratorSize(::Type{<:LazyKuhnSimplexes}) = Base.SizeUnknown()
Base.eltype(::Type{LazyKuhnSimplexes{N}}) where N = Simplex{N}

function Base.iterate(iter::LazyKuhnSimplexes{N}, state=nothing) where N
    perm_iter = KuhnPermutationIterator(N)

    # Zustand: nichts = Start, sonst vorheriger Permutations-Zustand
    perm_state = state

    # Explizite Unterscheidung: Start oder Fortsetzung
    next = perm_state === nothing ? iterate(perm_iter) : iterate(perm_iter, perm_state)
    next === nothing && return nothing

    perm, new_perm_state = next

    indices = generate_kuhn_indices(iter.cell_origin_idx, perm)

    vertex_data = [get_vertex!(iter.cache, idx) for idx in indices]
    grads = [d[2] for d in vertex_data]

    hull = ConvexHullArray([Singleton(g) for g in grads])
    zero_vec = zero(SVector{N, Float64})

    if !(zero_vec ∈ hull)
        # Gepruned → rekursiv mit nächstem Zustand (darf auch nothing sein)
        return iterate(iter, new_perm_state)
    end

    vertices = [
        iter.cache.lb .+ (SVector{N}(idx.I .- 1) .* iter.cache.cell_width)
        for idx in indices
    ]

    simplex = Simplex{N}(vertices, indices)

    return simplex, new_perm_state
end