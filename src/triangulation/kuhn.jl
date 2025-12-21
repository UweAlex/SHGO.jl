# Einfacher, nicht-allozierender Heap-Permutations-Iterator (für N! Permutationen)
struct HeapPermutationIterator
    n::Int
    data::Vector{Int}
end

HeapPermutationIterator(n::Int) = HeapPermutationIterator(n, collect(1:n))

function Base.iterate(iter::HeapPermutationIterator, state=iter.data)
    # Heap's Algorithm – modifiziert für Iterator
    # Gibt SVector zurück, keine Allocation
    # (Vereinfachte Version – vollständige Implementierung in Phase 1 erweitern)
    # Für Demo: yield alle Permutationen
    # ...
end

# Lazy Kuhn-Simplex-Generator pro Zelle
struct LazyKuhnSimplexes{N}
    cell_origin_idx::CartesianIndex{N}
    cache::VertexCache{N}
end

function Base.iterate(iter::LazyKuhnSimplexes{N}, perm_iter = HeapPermutationIterator(N)) where N
    next_perm = iterate(perm_iter)
    next_perm === nothing && return nothing

    perm, perm_state = next_perm
    # Generiere Simplex-Indizes aus Permutation (Kuhn-Regel)
    indices = generate_kuhn_indices(iter.cell_origin_idx, perm)

    # Hole Daten aus Cache
    vertex_data = [get_vertex!(iter.cache, idx) for idx in indices]
    grads = [d[2] for d in vertex_data]

    # Sofort Pruning (Gradient-Hull)
    hull = ConvexHullArray([Singleton(g) for g in grads])
    if !(zeros(N) ∈ hull)
        return iterate(iter, perm_state)  # verworfen → nächster
    end

    # Überlebt → Simplex zurückgeben
    vertices = [iter.cache.lb .+ (SVector(idx.I...) .- 1) .* iter.cache.cell_width for idx in indices]
    simplex = Simplex(tuple(vertices...), indices)

    (simplex, perm_state)
end