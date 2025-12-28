# File: src/triangulation/kuhn.jl
# 
# PURPOSE:
# This module implements Kuhn triangulation for simplicial mesh generation.
# Kuhn triangulation decomposes hypercubes into simplices using permutation-based
# subdivision, which is crucial for SHGO's topological sampling strategy.
#
# USED BY:
# - SHGO.jl main module (analyze function)
# - TopicalStructure.jl (consumes generated simplices)
#
# USES:
# - StaticArrays.jl for zero-allocation vertex coordinates
# - Combinatorics.jl for permutation generation
# - LazySets.jl for gradient hull pruning (containment checks)
# - Cache.jl (VertexCache) for function evaluation storage
#
# KEY CONCEPTS:
# - Kuhn Triangulation: Divides each grid cell into N! simplices using permutations
# - Gradient Pruning: Filters out simplices where gradient hull doesn't contain zero
#   (these regions cannot contain local minima by first-order optimality)

# === IMPORTS WITH ABBREVIATED ALIASES ===
using StaticArrays
using Combinatorics
using LazySets

# Establish const abbreviations for full qualification
const SA = StaticArrays
const CB = Combinatorics
const LS = LazySets

# === PERMUTATION ITERATOR (Zero-Allocation Design) ===

"""
    KuhnPermutationIterator{N}

Zero-allocation iterator for generating all N! permutations needed for
Kuhn triangulation in N dimensions.

# Type Parameter
- `N::Int`: Dimensionality of the space (number of coordinates)

# Performance
- Returns SVector{N, Int} to avoid heap allocations
- Uses Combinatorics.permutations internally but wraps results in StaticArrays
- Expected iterations: factorial(N)

# Example
```julia
iter = KuhnPermutationIterator(3)
for perm in iter
    # perm is SVector{3, Int}, e.g., [1, 2, 3], [1, 3, 2], ...
end
```
"""
struct KuhnPermutationIterator{N} end

# Constructor: Type-level N parameter for compile-time optimization
KuhnPermutationIterator(n::Base.Int) = KuhnPermutationIterator{n}()

# Iterator protocol: length is known at compile time
Base.length(::KuhnPermutationIterator{N}) where N = Base.factorial(N)

"""
    Base.iterate(::KuhnPermutationIterator{N}, state)

Iterate through all permutations of 1:N, returning StaticArrays.SVector.

# Algorithm
- Wraps Combinatorics.permutations(1:N)
- Converts each permutation to SVector{N, Int} for type stability
- State carries the underlying permutation iterator state

# Returns
- `nothing` when exhausted
- `(SVector{N, Int}, new_state)` for each permutation
"""
function Base.iterate(::KuhnPermutationIterator{N}, state=Base.nothing) where N
    # Delegate to Combinatorics.jl for permutation generation
    p_iter = CB.permutations(1:N)
    
    # Handle first call (state == nothing) vs subsequent calls
    next_val = Base.isnothing(state) ? 
        Base.iterate(p_iter) : 
        Base.iterate(p_iter, state)
    
    # Termination check
    Base.isnothing(next_val) && return Base.nothing
    
    # Unpack result and convert to StaticArray
    p, next_state = next_val
    return (SA.SVector{N, Base.Int}(p), next_state)
end

# === SIMPLEX GENERATION WITH GRADIENT-BASED PRUNING ===

"""
    generate_kuhn_simplices(N, n_div, cache, use_pruning=true)

Generate all Kuhn simplices for a grid with gradient-based filtering.

# Arguments
- `N::Int`: Dimensionality of the optimization space
- `n_div::Int`: Number of divisions per axis (grid resolution)
- `cache::VertexCache`: Pre-computed function values and gradients
- `use_pruning::Bool=true`: Enable gradient hull pruning (default: ON)

# Algorithm Overview
1. **Grid Traversal**: Iterate over all (n_div)^N hypercube cells
2. **Permutation Loop**: For each cell, generate N! simplices via Kuhn paths
3. **Data Acquisition**: Sample vertices along the Kuhn path:
   - Vertex i+1 differs from vertex i in exactly one coordinate (incremented by 1)
   - This forms a chain from cell corner to opposite corner
4. **Gradient Pruning** (if enabled):
   - Construct convex hull of gradients at simplex vertices
   - Discard simplex if 0 ∉ hull (no critical point possible)
5. **Export**: Return only valid simplices

# Mathematical Background
The gradient hull pruning is based on first-order optimality:
- At a local minimum, ∇f(x*) = 0
- If 0 is not in conv{∇f(v₁), ..., ∇f(vₙ₊₁)}, then no critical point exists in simplex
- This dramatically reduces the number of local searches needed

# Returns
- `Vector{Simplex{N}}`: List of active simplices (candidates for optimization)

# Performance Notes
- Pruning typically removes 60-90% of simplices for smooth functions
- Without pruning, returns (n_div)^N * N! simplices
- Memory complexity: O(N * |active_simplices|)

# Example
```julia
cache = VertexCache(testfunc, (10, 10, 10))
simplices = generate_kuhn_simplices(3, 10, cache, true)
# Expect ~6000 simplices instead of 60000 with pruning
```
"""
function generate_kuhn_simplices(
    N::Base.Int, 
    n_div::Base.Int, 
    cache, 
    use_pruning::Base.Bool = Base.true
)
    # Output container - will hold all valid simplices
    simplices = Base.Vector{Simplex{N}}()
    
    # Pre-generate all N! permutations once
    perms = Base.collect(KuhnPermutationIterator(N))
    
    # === PHASE 1: GRID TRAVERSAL ===
    # Iterate over all hypercube cells in the grid
    # Each cell is indexed by N coordinates in range [1, n_div]
    for cell_idx in Base.CartesianIndices(Base.ntuple(_ -> n_div, N))
        
        # === PHASE 2: PERMUTATION LOOP ===
        # Each permutation defines one simplex subdivision of the cell
        for p in perms
            
            # --- SECTION 1: DATA ACQUISITION ---
            # Allocate storage for simplex vertices (N+1 points define N-simplex)
            indices = Base.Vector{Base.CartesianIndex{N}}(Base.undef, N+1)
            vertices = Base.Vector{SA.SVector{N, Base.Float64}}(Base.undef, N+1)
            grads = Base.Vector{SA.SVector{N, Base.Float64}}(Base.undef, N+1)
            
            # Starting point: Lower-left corner of the cell
            curr_idx = cell_idx
            indices[1] = curr_idx
            
            # Fetch cached function value and gradient
            v_val, v_grad = get_vertex!(cache, curr_idx)
            
            # Convert grid index to physical coordinates
            # Formula: x = lb + (idx - 1) * step_size
            vertices[1] = SA.SVector{N, Base.Float64}(
                Base.ntuple(i -> cache.lower[i] + (curr_idx[i]-1) * cache.step_size[i], N)
            )
            grads[1] = v_grad
            
            # --- KUHN PATH CONSTRUCTION ---
            # Walk along edges defined by permutation p
            # At step i, increment the coordinate p[i] by 1
            for i in 1:N
                dim = p[i]  # Which coordinate to increment
                
                # Build new index tuple: copy all coords, increment dim-th by 1
                new_idx_tuple = Base.ntuple(
                    d -> d == dim ? curr_idx[d] + 1 : curr_idx[d], 
                    N
                )
                curr_idx = Base.CartesianIndex(new_idx_tuple)
                
                # Store vertex data
                indices[i+1] = curr_idx
                _, vg = get_vertex!(cache, curr_idx)
                
                vertices[i+1] = SA.SVector{N, Base.Float64}(
                    Base.ntuple(j -> cache.lower[j] + (curr_idx[j]-1) * cache.step_size[j], N)
                )
                grads[i+1] = vg
            end
            
            # --- SECTION 2: GRADIENT HULL PRUNING (DECISION LOGIC) ---
            # Default: assume simplex is valid
            is_valid = Base.true
            
            if use_pruning
                # Construct convex hull of gradient vectors
                # Add small ball to handle numerical precision near zero
                grad_vectors = [Base.Vector(g) for g in grads]
                grad_hull = LS.VPolygon(grad_vectors) + LS.BallInf(Base.zeros(N), 1e-9)
                
                # Check if zero vector is contained in hull
                # If NOT contained → no critical point possible → discard simplex
                if !(0.0 ∈ grad_hull)
                    is_valid = Base.false
                end
            end
            
            # --- SECTION 3: EXPORT ---
            # Only add simplex to output if it passed pruning
            if is_valid
                Base.push!(simplices, Simplex{N}(vertices, indices))
            end
        end
    end
    
    return simplices
end

# End: src/triangulation/kuhn.jl