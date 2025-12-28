# File: src/types.jl
#
# PURPOSE:
# Core type definitions for SHGO.jl geometric and result structures.
# Provides immutable data containers for simplicial complexes, spatial regions,
# and optimization results with solver-agnostic accessor functions.
#
# USED BY:
# - SHGO.jl main module (analyze function returns SHGOResult)
# - kuhn.jl (constructs Simplex objects during triangulation)
# - Basins.jl (clusters Simplex objects into Regions)
# - User-facing code (queries results via accessor functions)
#
# USES:
# - StaticArrays.jl for memory-efficient vertex coordinates
# - Base standard library for core types (Vector, CartesianIndex)
# - No runtime dependencies on optimization solvers (type-agnostic design)
#
# DESIGN PHILOSOPHY:
# - Immutable structs for thread-safety and cache optimization
# - Type parameters (N) for compile-time dimension specialization
# - Solver-agnostic result storage (Any type) for cross-version compatibility
# - Accessor pattern decouples user code from solver implementation details

# === IMPORTS WITH ABBREVIATED ALIASES ===
using StaticArrays
using NonlinearOptimizationTestFunctions
using SciMLBase

# Establish const abbreviations for full qualification
const SA = StaticArrays
const NOTF = NonlinearOptimizationTestFunctions
const SCIML = SciMLBase

# === GEOMETRIC PRIMITIVES ===

"""
    Simplex{N}

Immutable representation of an N-dimensional simplex (N+1 vertices).

# Type Parameter
- `N::Int`: Spatial dimensionality (e.g., N=2 for triangles, N=3 for tetrahedra)

# Fields
- `vertices::Vector{SVector{N, Float64}}`: Physical coordinates of vertices in ℝᴺ
- `indices::Vector{CartesianIndex{N}}`: Grid indices for cache lookups

# Mathematical Definition
A simplex is the convex hull of N+1 affinely independent points.
In 2D: triangle (3 vertices)
In 3D: tetrahedron (4 vertices)
In ND: N+1 vertices defining the minimal convex polytope

# Storage Strategy
- **vertices**: Actual coordinates for geometric queries (distance, containment)
- **indices**: Grid references for efficient value/gradient cache access
- Both stored redundantly to optimize different query patterns

# Invariants
- `length(vertices) == length(indices) == N+1`
- All vertices must be distinct (non-degenerate simplex)
- Indices reference valid grid cells within domain bounds

# Example
```julia
# 2D triangle simplex
s = Simplex{2}(
    [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(0.0, 1.0)],
    [CartesianIndex(1,1), CartesianIndex(2,1), CartesianIndex(1,2)]
)
```
"""
struct Simplex{N}
    vertices::Base.Vector{SA.SVector{N, Base.Float64}}
    indices::Base.Vector{Base.CartesianIndex{N}}
end

"""
    Region{N}

Immutable collection of simplices forming a connected topological region.

# Type Parameter
- `N::Int`: Spatial dimensionality (inherited from constituent simplices)

# Fields
- `simplices::Vector{Simplex{N}}`: Set of simplices comprising this region

# Topological Interpretation
A Region represents a "basin of attraction" or "star domain" in SHGO terminology:
- All simplices share topological connectivity (adjacent faces)
- All points in region's simplices converge to same local minimum under gradient flow
- Used for clustering simplices before local optimization

# Usage in SHGO Pipeline
1. **Triangulation** generates individual simplices
2. **Clustering** groups simplices into Regions based on:
   - Topological adjacency (shared facets)
   - Energy barrier analysis (gradient continuity)
3. **Optimization** launches one local search per Region

# Example
```julia
# Region containing 3 adjacent triangles
r = Region{2}([simplex1, simplex2, simplex3])
```
"""
struct Region{N}
    simplices::Base.Vector{Simplex{N}}
end

# === OPTIMIZATION RESULT CONTAINER ===

"""
    SHGOResult

Immutable container for SHGO analysis results with solver-agnostic storage.

# Design Rationale
The fields `global_minimum` and `local_minima` use `Any` typing to decouple
SHGO from specific solver implementations. Different optimization backends
(OptimizationOptimJL, NLopt, etc.) return different solution types, and
these types may change across package versions. By using `Any`:
- Cross-version compatibility is maintained
- Multiple solvers can be swapped without code changes
- Duck typing allows uniform access via `.objective` and `.u` properties

# Fields
- `global_minimum::Any`: Best solution found (or `nothing` if search failed)
- `local_minima::Vector{Any}`: All distinct local minima discovered
- `num_basins::Int`: Count of topologically distinct attraction basins

# Accessor Pattern (RECOMMENDED)
Do NOT access fields directly. Use provided accessor functions:
- `get_global_value(result)` → objective value
- `get_global_point(result)` → minimizer coordinates
- `get_local_values(result)` → all objective values
- `get_local_points(result)` → all minimizer coordinates

This pattern future-proofs code against solver API changes.

# Invariants
- `num_basins >= length(local_minima)` (some basins may share same minimum)
- If `global_minimum !== nothing`, it equals `argmin(local_minima)`
- All elements in `local_minima` satisfy domain bounds

# Example
```julia
result = analyze(test_function, n_div=20)
println("Global minimum: ", get_global_value(result))
println("Found ", result.num_basins, " basins")
```
"""
struct SHGOResult
    "Best solution found (Optimization.solve result or nothing)"
    global_minimum::Base.Any

    "All distinct local minima (Vector of solver solution objects)"
    local_minima::Base.Vector{Base.Any}

    "Number of topologically distinct basins of attraction"
    num_basins::Base.Int
end

# ═══════════════════════════════════════════════════════════════
# ACCESSOR FUNCTIONS (Future-Proof Interface)
# ═══════════════════════════════════════════════════════════════
# These functions abstract away solver-specific details and provide
# a stable API regardless of underlying optimization backend changes.

"""
    get_global_value(result::SHGOResult) -> Float64

Extract objective function value of the global minimum.

# Returns
- `Inf` if no global minimum was found (`global_minimum === nothing`)
- `result.global_minimum.objective` otherwise

# Notes
- Assumes solver solutions implement `.objective` property (SciML standard)
- Returns `Inf` as sentinel value to avoid `nothing` propagation in comparisons

# Example
```julia
result = analyze(rosenbrock, n_div=15)
f_best = get_global_value(result)
@assert f_best ≈ 0.0 atol=1e-6  # Known Rosenbrock minimum
```
"""
get_global_value(r::SHGOResult) = 
    Base.isnothing(r.global_minimum) ? Base.Inf : r.global_minimum.objective

"""
    get_global_point(result::SHGOResult) -> Vector{Float64}

Extract minimizer coordinates of the global minimum.

# Returns
- Empty vector `Float64[]` if no global minimum found
- `result.global_minimum.u` otherwise (SciML-style minimizer field)

# Notes
- `.u` is the SciML standard field name for solution points
- Returns empty vector (not `nothing`) for type stability in downstream code

# Example
```julia
result = analyze(sphere, n_div=10)
x_star = get_global_point(result)
@assert norm(x_star) < 1e-6  # Sphere minimum at origin
```
"""
get_global_point(r::SHGOResult) = 
    Base.isnothing(r.global_minimum) ? Base.Float64[] : r.global_minimum.u

"""
    get_local_values(result::SHGOResult) -> Vector{Float64}

Extract objective values of all local minima.

# Returns
- Vector of function values at each local minimum
- Empty if `local_minima` is empty

# Notes
- Useful for landscape analysis (multimodality, ruggedness metrics)
- Values are **not** guaranteed to be sorted (use `sort!` if needed)

# Example
```julia
result = analyze(himmelblau, n_div=20)
f_vals = get_local_values(result)
@assert length(f_vals) == 4  # Himmelblau has 4 minima
```
"""
get_local_values(r::SHGOResult) = 
    [m.objective for m in r.local_minima]

"""
    get_local_points(result::SHGOResult) -> Vector{Vector{Float64}}

Extract minimizer coordinates of all local minima.

# Returns
- Vector of coordinate vectors (one per local minimum)
- Empty if `local_minima` is empty

# Notes
- Each element is a `Vector{Float64}` of length N (dimensionality)
- Points are **not** guaranteed to be sorted by objective value

# Example
```julia
result = analyze(sixhump, n_div=15)
points = get_local_points(result)
for (i, x) in enumerate(points)
    println("Minimum ", i, " at ", x)
end
```
"""
get_local_points(r::SHGOResult) = 
    [m.u for m in r.local_minima]

# === EXPORTS ===
# Public API: types and accessor functions
Base.export(
    Simplex, 
    Region, 
    SHGOResult,
    get_global_value, 
    get_global_point,
    get_local_values, 
    get_local_points
)

# End: src/types.jl