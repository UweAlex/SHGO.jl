# File: src/cache.jl
#
# PURPOSE:
# Thread-safe memoization layer for expensive function evaluations during SHGO sampling.
# Prevents redundant calls to objective/gradient functions by caching results indexed
# by grid coordinates. Critical for performance when the same vertices are accessed
# multiple times during triangulation and refinement phases.
#
# USED BY:
# - kuhn.jl (queries vertices during simplex generation)
# - SHGO.jl main module (initializes cache before triangulation)
# - TopicalStructure.jl (accesses cached gradients for pruning)
#
# USES:
# - StaticArrays.jl for memory-efficient coordinate storage
# - Base.Threads for concurrent access synchronization (ReentrantLock)
# - NonlinearOptimizationTestFunctions.jl for test function interface
#
# KEY FEATURES:
# - Thread-safe lazy evaluation (lock-protected Dict access)
# - Floating-point precision safeguards (direct coordinate calculation)
# - Bounds validation to catch indexing errors early
# - Optional duplicate detection for debugging triangulation issues
#
# PERFORMANCE NOTES:
# - Cache hit: O(1) dictionary lookup (no function evaluation)
# - Cache miss: O(f) where f = cost of objective + gradient evaluation
# - Typical hit rate: 60-80% for Kuhn triangulation (vertices shared between simplices)
# - Memory usage: O(N * |vertices|) where N = dimensionality

module Cache

# === IMPORTS WITH ABBREVIATED ALIASES ===
using StaticArrays
using Base.Threads: ReentrantLock
using NonlinearOptimizationTestFunctions: TestFunction

# Establish const abbreviations for full qualification
const SA = StaticArrays
const BT = Base.Threads
const NOTF = NonlinearOptimizationTestFunctions

# === PUBLIC API EXPORTS ===
Base.export(VertexCache, get_vertex!)

# === CACHE DATA STRUCTURE ===

"""
    VertexCache{N, T<:TestFunction}

Thread-safe memoization cache for vertex function evaluations on a grid.

# Type Parameters
- `N::Int`: Dimensionality of the optimization space
- `T<:TestFunction`: Concrete type of the test function being optimized

# Fields
- `storage::Dict{CartesianIndex{N}, Tuple{Float64, SVector{N, Float64}}}`: 
  Hash map from grid indices to (value, gradient) pairs
- `lock::ReentrantLock`: Synchronization primitive for concurrent access
- `tf::T`: Test function instance providing f(x) and ∇f(x)
- `lb::SVector{N, Float64}`: Lower bounds of search domain (per dimension)
- `ub::SVector{N, Float64}`: Upper bounds of search domain (per dimension)
- `cell_width::SVector{N, Float64}`: Grid spacing (precomputed: (ub - lb) / divisions)

# Design Decisions

## Why CartesianIndex as Key?
Grid indices are integers → no floating-point comparison issues.
Two vertices with "close" coordinates might have different indices due to
numerical errors, but two accesses to the **same grid cell** always hit cache.

## Why Store Both Value and Gradient?
SHGO needs both for:
1. Value-based pruning (eliminate high-value regions)
2. Gradient-based pruning (eliminate regions where ∇f ≠ 0 everywhere)
Fetching separately would require two function calls per vertex.

## Why ReentrantLock?
- Allows same thread to acquire lock multiple times (recursive calls safe)
- Standard Mutex would deadlock if get_vertex! calls itself
- Small overhead (~10ns) compared to function evaluation cost (~μs to ms)

# Thread Safety Guarantee
Multiple threads can safely call `get_vertex!` concurrently on the same cache.
The lock ensures only one thread evaluates a given vertex, others wait and
receive the cached result.

# Example
```julia
tf = TEST_FUNCTIONS["rosenbrock"]
cache = VertexCache(tf, (20, 20))  # 20x20 grid in 2D

# Thread-safe access from multiple tasks
@threads for i in 1:100
    idx = CartesianIndex(rand(1:20), rand(1:20))
    val, grad = get_vertex!(cache, idx)
end
```
"""
struct VertexCache{N, T<:NOTF.TestFunction}
    storage::Base.Dict{Base.CartesianIndex{N}, Base.Tuple{Base.Float64, SA.SVector{N, Base.Float64}}}
    lock::BT.ReentrantLock
    tf::T
    lb::SA.SVector{N, Base.Float64}
    ub::SA.SVector{N, Base.Float64}
    cell_width::SA.SVector{N, Base.Float64}
end

# === CONSTRUCTOR WITH BOUNDS OVERRIDE ===

"""
    VertexCache(tf::TestFunction, divisions::NTuple{N, Int}; lb=nothing, ub=nothing) where N

Construct a vertex cache for a test function with optional custom bounds.

# Arguments
- `tf::TestFunction`: Test function to be optimized (must implement `.f` and `.grad`)
- `divisions::NTuple{N, Int}`: Number of grid cells per dimension (not vertices!)
- `lb::Union{Nothing, AbstractVector}`: Custom lower bounds (default: use `tf.lb`)
- `ub::Union{Nothing, AbstractVector}`: Custom upper bounds (default: use `tf.ub`)

# Algorithm
1. Resolve bounds: use custom if provided, otherwise extract from test function
2. Compute cell width: `(ub - lb) / divisions` per dimension
3. Allocate empty storage dictionary (lazy population during access)
4. Initialize reentrant lock for thread safety

# Why Allow Bounds Override?
- Testing: Create small caches with known coordinate ranges
- Domain restriction: Focus search on subset of function's natural domain
- Multi-scale analysis: Compare behavior on different scales

# Bounds Validation
Constructor does NOT validate that lb < ub (deferred to first get_vertex! call).
This allows degenerate caches for testing edge cases.

# Example
```julia
# Standard construction (use function's natural bounds)
tf = fixed(TEST_FUNCTIONS["sphere"], n=3)
cache1 = VertexCache(tf, (10, 10, 10))

# Custom bounds (zoom into region near origin)
cache2 = VertexCache(tf, (50, 50, 50), lb=[-0.1, -0.1, -0.1], ub=[0.1, 0.1, 0.1])
```
"""
function VertexCache(
    tf::NOTF.TestFunction, 
    divisions::Base.NTuple{N, Base.Int}; 
    lb=Base.nothing, 
    ub=Base.nothing
) where N
    # Resolve bounds with fallback to test function defaults
    actual_lb = Base.isnothing(lb) ? 
        SA.SVector{N}(tf.lb) : 
        SA.SVector{N}(lb)
    
    actual_ub = Base.isnothing(ub) ? 
        SA.SVector{N}(tf.ub) : 
        SA.SVector{N}(ub)
    
    # Precompute cell width for efficient coordinate calculation
    # divisions = number of cells, not vertices (vertices = divisions + 1)
    width = (actual_ub .- actual_lb) ./ SA.SVector{N}(divisions)
    
    # Construct cache with empty storage (lazy evaluation pattern)
    return VertexCache(
        Base.Dict{Base.CartesianIndex{N}, Base.Tuple{Base.Float64, SA.SVector{N, Base.Float64}}}(),
        BT.ReentrantLock(),
        tf,
        actual_lb,
        actual_ub,
        width
    )
end

# === THREAD-SAFE VERTEX ACCESS ===

"""
    get_vertex!(cache::VertexCache{N}, idx::CartesianIndex{N}; warn_duplicates=true) -> (Float64, SVector{N, Float64})

Retrieve or compute function value and gradient at grid vertex.

# Arguments
- `cache::VertexCache{N}`: Cache instance to query
- `idx::CartesianIndex{N}`: Grid coordinates (1-indexed Julia convention)
- `warn_duplicates::Bool=true`: Emit warning on cache hit (for debugging)

# Returns
- `(value::Float64, gradient::SVector{N, Float64})`: Function evaluation results

# Algorithm (Thread-Safe Lazy Evaluation)
1. **Acquire lock**: Block if another thread is accessing storage
2. **Check cache**: If `idx` exists in storage, return immediately (cache hit)
3. **Compute coordinates**: `pos = lb + (idx - 1) * cell_width` (vectorized)
4. **Validate bounds**: Error if computed position exceeds domain limits
5. **Evaluate function**: Call `tf.f(pos)` and `tf.grad(pos)`
6. **Store result**: Insert into cache for future hits
7. **Release lock**: Allow other threads to proceed

# Floating-Point Precision Strategy
Uses **direct formula** instead of cumulative addition to avoid error accumulation:
- ❌ Bad: `pos = lb; for i in 1:n: pos += step` → O(n) rounding errors
- ✅ Good: `pos = lb + (idx-1) * step` → O(1) rounding error

# Bounds Validation
Catches two classes of errors:
1. **Index out of range**: User provided idx > grid dimensions
2. **Numerical overflow**: Floating-point arithmetic produced pos > ub

Both indicate bugs in triangulation code and should fail fast.

# Thread Safety Proof
- Lock acquired before any shared state access (storage Dict)
- `get!` is atomic within lock scope (no race between check and insert)
- Lock released automatically on block exit (even if exception thrown)

# Performance Characteristics
- **Cache hit**: O(1) hash lookup + O(1) lock acquisition ≈ 50-100ns
- **Cache miss**: O(f) where f = cost of objective + gradient ≈ 1μs-1ms
- **Hit rate**: Typically 60-80% for Kuhn triangulation

# Example
```julia
cache = VertexCache(test_func, (100, 100))

# First access: computes and stores
val1, grad1 = get_vertex!(cache, CartesianIndex(50, 50))  # ~1μs (miss)

# Second access: instant retrieval
val2, grad2 = get_vertex!(cache, CartesianIndex(50, 50))  # ~50ns (hit)

@assert val1 == val2  # Deterministic caching
```
"""
function get_vertex!(
    cache::VertexCache{N}, 
    idx::Base.CartesianIndex{N}; 
    warn_duplicates::Base.Bool=Base.true
) where N
    # Thread-safe critical section
    Base.lock(cache.lock) do
        # Duplicate detection (optional debugging feature)
        if Base.haskey(cache.storage, idx) && warn_duplicates
            Base.@warn "Duplicate index detected: $idx – Deduplication active (skipping computation)."
        end
        
        # Lazy evaluation: compute only if not cached
        Base.get!(cache.storage, idx) do
            # === COORDINATE CALCULATION (Precision-Critical) ===
            # Convert CartesianIndex to 0-indexed tuple for arithmetic
            idx_tuple = idx.I .- 1
            
            # Direct formula avoids cumulative floating-point errors
            pos = cache.lb .+ idx_tuple .* cache.cell_width
            
            # === BOUNDS VALIDATION ===
            # Catch indexing bugs early (before expensive function call)
            if Base.any(pos .> cache.ub) || Base.any(pos .< cache.lb)
                Base.error(
                    "Index $idx exceeds bounds: pos=$pos " *
                    "(lb=$(cache.lb), ub=$(cache.ub))"
                )
            end
            
            # === FUNCTION EVALUATION ===
            # This is the expensive operation being cached
            val = cache.tf.f(pos)
            grad = cache.tf.grad(pos)
            
            # Return tuple for storage
            (val, grad)
        end
    end
end

# === TESTING UTILITIES ===

"""
    force_duplicate_for_test!(cache::VertexCache{N}, idx::CartesianIndex{N})

Inject fake cache entry to test duplicate detection logic (testing only).

# Purpose
Allows unit tests to verify that:
1. Duplicate warnings are emitted correctly
2. Cached values are returned (not recomputed)
3. Race conditions don't corrupt storage

# Warning
**Never use in production code.** This function intentionally corrupts cache
with dummy values (0.0 everywhere) to trigger duplicate detection paths.

# Example
```julia
@testset "Duplicate Detection" begin
    cache = VertexCache(test_func, (10, 10))
    idx = CartesianIndex(5, 5)
    
    # Force fake cache entry
    force_duplicate_for_test!(cache, idx)
    
    # Next access should warn
    @test_logs (:warn, r"Duplicate index") get_vertex!(cache, idx)
end
```
"""
function force_duplicate_for_test!(
    cache::VertexCache{N}, 
    idx::Base.CartesianIndex{N}
) where N
    # Compute position (for realistic dummy entry)
    pos = cache.lb .+ (idx.I .- 1) .* cache.cell_width
    
    # Manually inject dummy entry
    cache.storage[idx] = (0.0, SA.SVector{N}(0.0))
    
    # Trigger duplicate warning path
    get_vertex!(cache, idx)
end

end # module Cache

# End: src/cache.jl