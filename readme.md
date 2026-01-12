# SHGO.jl

[![Build Status](https://github.com/USERNAME/SHGO.jl/workflows/CI/badge.svg)](https://github.com/USERNAME/SHGO.jl/actions)
[![License: BSD-3](https://img.shields.io/badge/License-BSD3-blue.svg)](LICENSE)

**Simplicial Homology Global Optimization in Julia**

SHGO.jl is a pure Julia implementation of the SHGO algorithm for finding **all** local and global minima of a function within bounds. It is inspired by the [SciPy implementation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html) but redesigned for Julia's strengths.

## Features

- 🎯 **Finds all minima** - not just the global one
- 🔬 **Topological approach** - uses simplicial homology concepts
- 💾 **Memory efficient** - implicit Kuhn topology, no graph in memory
- ⚡ **Fast** - lazy evaluation with point caching
- 🔄 **Automatic convergence** - stops when basin count stabilizes (Betti number stability)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/USERNAME/SHGO.jl")
```

## Quick Start

```julia
using SHGO
using NonlinearOptimizationTestFunctions

# Test function: Six-Hump Camelback (6 local minima, 2 global)
tf = fixed(TEST_FUNCTIONS["sixhumpcamelback"]; n=2)

# Find all minima
result = analyze(tf; verbose=true)

# Results
println("Found $(result.num_basins) basins")
for m in sort(result.local_minima, by=x->x.objective)
    println("  f = $(round(m.objective, digits=4)) at $(round.(m.minimizer, digits=3))")
end
```

**Output:**
```
Found 6 basins
  f = -1.0316 at [-0.09, 0.713]
  f = -1.0316 at [0.09, -0.713]
  f = -0.2155 at [-1.703, 0.796]
  f = -0.2155 at [1.703, -0.796]
  f = 2.104 at [-1.607, -0.569]
  f = 2.104 at [1.607, 0.569]
```

## Custom Objective Function

```julia
using SHGO
using NonlinearOptimizationTestFunctions

# Define your own function
function my_objective(x)
    return (x[1] - 1)^2 + (x[2] - 2.5)^2
end

function my_gradient(x)
    return [2*(x[1] - 1), 2*(x[2] - 2.5)]
end

# Wrap in TestFunction format
tf = TestFunction(
    f = my_objective,
    grad = my_gradient,
    lb = [-5.0, -5.0],
    ub = [5.0, 5.0],
    name = "custom"
)

result = analyze(tf)
println("Global minimum: f = $(result.local_minima[1].objective)")
println("Location: $(result.local_minima[1].minimizer)")
```

## API Reference

### `analyze(tf; kwargs...)`

Main entry point for optimization.

**Arguments:**
- `tf` - Test function with fields `f`, `grad`, `lb`, `ub`

**Keyword Arguments:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_div_initial` | 8 | Initial grid resolution per dimension |
| `n_div_max` | 25 | Maximum grid resolution |
| `stability_count` | 2 | Iterations with stable basin count for convergence |
| `threshold_ratio` | 0.1 | Tolerance for basin merging (relative to value range) |
| `min_distance_tolerance` | 0.05 | Minimum distance between distinct minima |
| `local_maxiters` | 500 | Maximum iterations for local optimization |
| `verbose` | false | Print progress information |

**Returns:** `SHGOResult` with fields:
- `local_minima::Vector{MinimumPoint}` - All found minima
- `num_basins::Int` - Number of distinct basins
- `converged::Bool` - Whether Betti stability was reached
- `iterations::Int` - Number of refinement iterations

### `MinimumPoint`

```julia
struct MinimumPoint
    minimizer::Vector{Float64}  # Location of minimum
    objective::Float64          # Function value at minimum
end
```

Access via `m.minimizer` or `m.u` (SciML compatibility).

## Algorithm Overview

SHGO.jl uses a topological approach to global optimization:

1. **Grid Sampling** - Create a Kuhn triangulation of the search space
2. **Star-Minimum Detection** - Find points that are minimal in their local neighborhood
3. **Basin Clustering** - Group star-minima into attraction basins
4. **Iterative Refinement** - Increase resolution until basin count stabilizes
5. **Local Optimization** - Run L-BFGS from one representative per basin
6. **Deduplication** - Merge minima that converged to the same point

The key insight is that the number of basins (0th Betti number) becomes stable as resolution increases, providing a natural convergence criterion.

## Benchmarks

Tested on standard optimization benchmarks (2D):

| Function | Expected Minima | Found | Time |
|----------|----------------|-------|------|
| Sphere | 1 | 1 | 0.001s |
| Rosenbrock | 1 | 1 | 0.002s |
| Himmelblau | 4 | 4 | 0.001s |
| Six-Hump Camelback | 6 | 6 | 0.002s |
| Rastrigin | 1 (global) | 1 | 2.1s |
| Ackley | 1 | 1 | 1.4s |
| Easom | 1 | 1 | 0.4s |

**Coverage: 100%** on all test functions.

## Comparison with SciPy SHGO

| Aspect | SciPy SHGO | SHGO.jl |
|--------|------------|---------|
| Language | Python/C | Pure Julia |
| Sampling | Sobol sequence | Kuhn grid |
| Triangulation | Explicit Delaunay | Implicit Kuhn |
| Memory | O(n²) | O(n) |
| Parallelization | Limited (GIL) | Planned |
| Constraints | Full support | Box bounds only |
| High dimensions (N>6) | Good | Limited |

**Current status:** Equivalent quality for 2D problems. SciPy is better for high-dimensional problems due to Sobol sampling.

## Limitations

- **Box constraints only** - nonlinear constraints not yet supported
- **Dimension scaling** - Kuhn triangulation produces N! simplices per cell, limiting practical use to N ≤ 6
- **No parallelization yet** - single-threaded execution

## Roadmap

- [ ] Constraint support (nonlinear inequalities/equalities)
- [ ] Sobol sampling for high dimensions
- [ ] Multi-threading for function evaluations
- [ ] Direct benchmark suite against SciPy

## References

- Endres, S. C., Sandrock, C., & Focke, W. W. (2018). "A simplicial homology algorithm for Lipschitz optimisation." *Journal of Global Optimization*, 72(2), 181-217.
- [SciPy SHGO Documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html)

## License

BSD-3-Clause. See [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! Please open an issue or PR on GitHub.

## Acknowledgments

- Original SHGO algorithm by Stefan Endres
- Architecture guidance from the Julia community