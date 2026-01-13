# bench.jl – 30 Funktionen, 1 Run, mit f/g-Calls & Coverage %
# Stand: 02.01.2026 – NUR bukin5 durch michalewicz ersetzt, sonst 100 % deine Version
import Pkg; Pkg.add("PyCall")
using SHGO
using NonlinearOptimizationTestFunctions
using PyCall
using Printf
using BenchmarkTools
using Logging
using LinearAlgebra

global_logger(SimpleLogger(stderr, Logging.Error))

const NOTF = NonlinearOptimizationTestFunctions
const SP   = pyimport("scipy.optimize")

# 30 klassische Testfunktionen (2D) – alle sofort nutzbar
const FUNCTIONS = [
    "sphere", "rosenbrock", "ackley", "rastrigin", "griewank", "levy",
    "beale", "himmelblau", "goldsteinprice", "bukin6", "threehumpcamel",
    "easom", "crossintray", "eggholder", "dropwave", "holdertable",
    "mccormick", "schwefel", "booth", "matyas", "bukin2",
    "bird", "brent", "michalewicz", "xinsheyang2", "xinsheyang3",
    "xinsheyang4", "schaffer2", "schaffer4", "sixhumpcamelback"
]

# Erweiterte bekannte Minima (global + lokal, aus Literatur)
const KNOWN_MINIMA = Dict(
    "sphere"           => [([0.0, 0.0], 0.0)],
    "rosenbrock"       => [([1.0, 1.0], 0.0)],
    "ackley"           => [([0.0, 0.0], 0.0)],
    "rastrigin"        => [([0.0, 0.0], 0.0)],
    "griewank"         => [([0.0, 0.0], 0.0)],
    "levy"             => [([1.0, 1.0], 0.0)],
    "beale"            => [([3.0, 0.5], 0.0)],
    "himmelblau"       => [
        ([3.0, 2.0], 0.0), ([-2.805, 3.131], 0.0),
        ([-3.779, -3.283], 0.0), ([3.584, -1.848], 0.0)
    ],
    "goldsteinprice"   => [([0.0, -1.0], 3.0)],
    "bukin6"           => [([-10.0, 1.0], 0.0)],
    "threehumpcamel"   => [([0.0, 0.0], 0.0)],
    "easom"            => [([π, π], -1.0)],
    "crossintray"      => [([1.3491, 1.3491], -2.0626)],  # eines von 4
    "eggholder"        => [([512.0, 404.2319], -959.6407)],
    "dropwave"         => [([0.0, 0.0], -1.0)],
    "holdertable"      => [([1.306, 1.306], -19.2085)],
    "mccormick"        => [([-0.54719, -1.54719], -1.9133)],
    "schwefel"         => [([420.9687, 420.9687], 0.0)],
    "booth"            => [([1.0, 3.0], 0.0)],
    "matyas"           => [([0.0, 0.0], 0.0)],
    "bukin2"           => [([-10.0, 1.0], 0.0)],
    "bird"             => [([4.701, 3.153], -106.765)],
    "brent"            => [([-10.0, 0.0], 0.0)],
    "michalewicz"      => [([2.20290552, 1.57079633], -1.8013)],  # Multimodal, global at (2.20, 1.57)
    "xinsheyang2"      => [([0.0, 0.0], 0.0)],
    "xinsheyang3"      => [([0.0, 0.0], -1.0)],
    "xinsheyang4"      => [([0.0, 0.0], -1.0)],
    "schaffer2"        => [([0.0, 0.0], 0.0)],
    "schaffer4"        => [([0.0, 0.0], 0.292579)],
    "sixhumpcamelback" => [
        ([-0.0898, 0.7126], -1.0316),
        ([0.0898, -0.7126], -1.0316),
        ([-1.7036, 0.7961], -0.2155),
        ([1.7036, -0.7961], -0.2155),
        ([-1.6071, -0.5687], 2.1040),
        ([1.6071, 0.5687], 2.1040)
    ]
)

function run_jl(tf, n_div; use_pruning=false, verbose=false)
    NOTF.reset_counts!(tf)

    res = with_logger(SimpleLogger(devnull, Logging.Error)) do
        SHGO.analyze(
            tf;
            n_div = n_div,
            use_gradient_pruning = use_pruning,
            verbose = verbose,
            adaptive_max_levels = 8,
            local_maxiters = 200
        )
    end

    return (
        result = res,
        f_calls = get_f_count(tf),
        g_calls = get_grad_count(tf)
    )
end

function run_py(tf, n_py)
    NOTF.reset_counts!(tf)

    f_py = x -> Float64(tf.f(x))
    g_py = x -> Vector{Float64}(tf.grad(x))

    bounds = [(Float64(NOTF.lb(tf)[i]), Float64(NOTF.ub(tf)[i]))
              for i in 1:length(NOTF.lb(tf))]

    res = SP.shgo(
        f_py,
        bounds;
        n = Int(n_py),
        iters = 3,
        sampling_method = "sobol",
        minimizer_kwargs = Dict("jac" => g_py)
    )

    return (
        result = res,
        f_calls = get_f_count(tf),
        g_calls = get_grad_count(tf)
    )
end

function quality(res, fname; pos_tol=0.05, is_py=false)
    if is_py
        xl = get(res, "xl", [])
        funl = get(res, "funl", [])

        if isempty(xl)
            local_minima = Vector{Float64}[]
            objectives = Float64[]
        else
            local_minima = [Vector{Float64}(row) for row in eachrow(xl)]
            objectives = [Float64(f) for f in funl]
        end

        best_value = Float64(get(res, "fun", Inf))
        minimizer = get(res, "x", nothing) isa Nothing ? nothing : Vector{Float64}(get(res, "x", []))
    else
        local_minima = [m.minimizer for m in res.local_minima]
        objectives = [m.objective for m in res.local_minima]
        best_value = isempty(objectives) ? Inf : minimum(objectives)
        minimizer = isempty(local_minima) ? nothing : local_minima[argmin(objectives)]
    end

    num_basins = length(objectives)

    # Echte Coverage: Wie viele bekannte Minima wurden gefunden?
    known_list = get(KNOWN_MINIMA, fname, [])
    found_count = 0

    for (known_pos, known_val) in known_list
        found = false
        for (i, pos) in enumerate(local_minima)
            if norm(pos - known_pos) < pos_tol && abs(objectives[i] - known_val) < 1e-3
                found = true
                break
            end
        end
        found_count += found ? 1 : 0
    end

    coverage_percent = length(known_list) > 0 ? (found_count / length(known_list)) * 100 : 0.0

    return (basins = num_basins,
            coverage_percent = coverage_percent)
end

function benchmark(n_div=10)
    n_py = (n_div + 1)^2

    println("="^200)
    println("SHGO Benchmark – 30 Funktionen, 1 Run, mit f/g-Calls & Coverage %")
    println("n_div = $n_div | SciPy n = $n_py")
    println("="^200)

    @printf(
        "%-25s | %10s | %10s | %12s | %12s | %12s | %12s | %12s | %12s | %7s\n",
        "Function", "Basins JL", "Basins PY", "Coverage JL %", "Coverage PY %",
        "f_calls JL", "g_calls JL", "f_calls PY", "g_calls PY", "Speedup"
    )
    println("-"^200)

    for fn in FUNCTIONS
        tf = NOTF.fixed(NOTF.TEST_FUNCTIONS[fn]; n=2)

        t_jl = @belapsed run_jl($tf, $n_div)
        jl = run_jl(tf, n_div)
        q_jl = quality(jl.result, fn)

        t_py = @belapsed run_py($tf, $n_py)
        py = run_py(tf, n_py)
        q_py = quality(py.result, fn; is_py=true)

        speedup = t_py / t_jl

        @printf(
            "%-25s | %10d | %10d | %12.1f%% | %12.1f%% | %12d | %12d | %12d | %12d | %6.2fx\n",
            fn, q_jl.basins, q_py.basins,
            q_jl.coverage_percent, q_py.coverage_percent,
            jl.f_calls, jl.g_calls,
            py.f_calls, py.g_calls,
            speedup
        )
    end

    println("="^200)
end

println("\nStarting 30-Function-Benchmark...\n")
benchmark()