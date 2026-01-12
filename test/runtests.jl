# File: test/runtests.jl
#
# Main test entry point for SHGO.jl
# Run with: julia --project=. -e 'using Pkg; Pkg.test()'
#

using Test
using SHGO
using NonlinearOptimizationTestFunctions
using LinearAlgebra

const NOTF = NonlinearOptimizationTestFunctions

# Known minima from literature
const KNOWN_MINIMA = Dict(
    "sphere" => [([0.0, 0.0], 0.0)],
    "rosenbrock" => [([1.0, 1.0], 0.0)],
    "himmelblau" => [
        ([3.0, 2.0], 0.0),
        ([-2.805118, 3.131312], 0.0),
        ([-3.779310, -3.283186], 0.0),
        ([3.584428, -1.848126], 0.0)
    ],
    "sixhumpcamelback" => [
        ([-0.0898, 0.7126], -1.0316),
        ([0.0898, -0.7126], -1.0316),
    ]
)

function count_found(result, known; pos_tol=0.15, val_tol=0.02)
    found = 0
    for (pos, val) in known
        for m in result.local_minima
            if norm(m.minimizer - pos) < pos_tol && abs(m.objective - val) < val_tol
                found += 1
                break
            end
        end
    end
    return found
end

@testset "SHGO.jl" begin

    @testset "Core Types" begin
        @test isdefined(SHGO, :MinimumPoint)
        @test isdefined(SHGO, :SHGOResult)
        @test isdefined(SHGO, :PointCache)
        @test isdefined(SHGO, :KuhnTopology)
    end

    @testset "PointCache" begin
        cache = SHGO.PointCache([0.0, 0.0], [1.0, 1.0], [10, 10])
        
        # Position calculation
        @test SHGO.index_to_position(cache, (0, 0)) ≈ [0.0, 0.0]
        @test SHGO.index_to_position(cache, (10, 10)) ≈ [1.0, 1.0]
        
        # Bounds check
        @test SHGO.is_valid_index(cache, (5, 5)) == true
        @test SHGO.is_valid_index(cache, (11, 5)) == false
        @test SHGO.is_valid_index(cache, (-1, 5)) == false
        
        # Caching behavior
        f = x -> sum(x.^2)
        val1 = SHGO.get_value!(cache, (5, 5), f)
        @test SHGO.num_evaluated(cache) == 1
        val2 = SHGO.get_value!(cache, (5, 5), f)
        @test SHGO.num_evaluated(cache) == 1  # No re-evaluation
        
        # Infinity padding
        @test SHGO.get_value!(cache, (100, 100), f) == Inf
    end

    @testset "KuhnTopology" begin
        topo = SHGO.KuhnTopology([10, 10])
        
        neighbors = SHGO.get_neighbors(topo, (5, 5))
        @test (4, 5) in neighbors
        @test (6, 5) in neighbors
        @test (5, 4) in neighbors
        @test (5, 6) in neighbors
        
        simplices = SHGO.get_simplices_in_cube(topo, (0, 0))
        @test length(simplices) == 2  # 2! = 2 in 2D
    end

    @testset "Single-Minimum Functions" begin
        for fname in ["sphere", "rosenbrock"]
            @testset "$fname" begin
                tf = NOTF.fixed(NOTF.TEST_FUNCTIONS[fname]; n=2)
                result = SHGO.analyze(tf; n_div_initial=10, n_div_max=20)
                
                @test result.num_basins >= 1
                @test !isempty(result.local_minima)
                
                known = KNOWN_MINIMA[fname]
                @test count_found(result, known) >= 1
            end
        end
    end

    @testset "Multi-Minimum Functions" begin
        @testset "Himmelblau (4 minima)" begin
            tf = NOTF.fixed(NOTF.TEST_FUNCTIONS["himmelblau"]; n=2)
            result = SHGO.analyze(tf; n_div_initial=12, n_div_max=25)
            
            @test result.num_basins >= 2
            known = KNOWN_MINIMA["himmelblau"]
            @test count_found(result, known) >= 2
        end

        @testset "Six-Hump Camelback (6 minima)" begin
            tf = NOTF.fixed(NOTF.TEST_FUNCTIONS["sixhumpcamelback"]; n=2)
            result = SHGO.analyze(tf; n_div_initial=12, n_div_max=30)
            
            @test result.num_basins >= 2
            
            # Must find global minimum
            best = minimum(m.objective for m in result.local_minima)
            @test best < -1.0
        end
    end

    @testset "Convergence" begin
        tf = NOTF.fixed(NOTF.TEST_FUNCTIONS["sphere"]; n=2)
        result = SHGO.analyze(tf; n_div_initial=5, n_div_max=20, stability_count=2)
        
        @test result.converged == true
        @test result.iterations <= 10
    end

    @testset "Backward Compatibility" begin
        tf = NOTF.fixed(NOTF.TEST_FUNCTIONS["sphere"]; n=2)
        
        # Old API should work
        result = SHGO.analyze(tf; n_div=10, use_gradient_pruning=false)
        
        @test result.num_basins >= 1
        @test !isempty(result.local_minima)
    end

end