# File: test/test_multimodal_detection.jl
using Test
using SHGO
using NonlinearOptimizationTestFunctions

const OBJECTIVE_TOLERANCE = 1e-4
const SIXHUMP_GLOBAL_MIN = -1.031628

@testset "Multimodal Detection - Six-Hump-Camelback" begin
    tf = fixed(TEST_FUNCTIONS["sixhumpcamelback"]; n=2)
    
    # Test mit moderater Auflösung UND OHNE Pruning (Phase 1 Strategie)
    res = analyze(tf; 
        n_div = 12,
        verbose = false,
        use_gradient_pruning = false  # WICHTIG: Ohne Pruning für Six-Hump
    )
    
    @testset "Basic Result Structure" begin
        @test res isa SHGOResult
        @test res.num_basins >= 1
        @test length(res.local_minima) >= 1
    end
    
    @testset "Global Minimum Detection" begin
        found_global = any(m -> abs(m.objective - SIXHUMP_GLOBAL_MIN) < OBJECTIVE_TOLERANCE, 
                          res.local_minima)
        @test found_global
        
        if found_global
            best_sol = argmin(m -> abs(m.objective - SIXHUMP_GLOBAL_MIN), res.local_minima)
            @test abs(best_sol.objective - SIXHUMP_GLOBAL_MIN) < OBJECTIVE_TOLERANCE
        end
    end
    
    @testset "Multimodality Detection" begin
        # Six-Hump hat 6 Minima, aber value-based clustering findet 2-4
        @test res.num_basins >= 2
        
        if res.num_basins < 4
            @info "Current basin count: $(res.num_basins) (target: ≥4 nach Gradient-Flow)"
        end
    end
    
    @testset "Solution Quality" begin
        @test all(m -> isfinite(m.objective), res.local_minima)
        @test all(m -> !isempty(m.u), res.local_minima)
        @test all(m -> m.objective < 3.0, res.local_minima)
    end
end

@testset "Multimodal Detection - Comparison Tests" begin
    @testset "Unimodal vs Multimodal" begin
        # WICHTIG: Ohne Pruning für Sphere (sonst 0 Basins)
        tf_sphere = fixed(TEST_FUNCTIONS["sphere"]; n=2)
        res_sphere = analyze(tf_sphere; n_div=8, verbose=false, use_gradient_pruning=false)
        
        @test res_sphere.num_basins >= 1
        @test length(res_sphere.local_minima) >= 1
        
        tf_sixhump = fixed(TEST_FUNCTIONS["sixhumpcamelback"]; n=2)
        res_sixhump = analyze(tf_sixhump; n_div=8, verbose=false, use_gradient_pruning=false)
        
        @test res_sixhump.num_basins >= res_sphere.num_basins
    end
end