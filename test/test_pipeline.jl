# File: test\test_pipeline.jl
using Test
using SHGO
using NonlinearOptimizationTestFunctions
using StaticArrays

@testset "Pipeline Tests" begin
    # 1. Test: Rosenbrock (Skalierbare Funktion aus der Library)
    # Wir holen sie sicher aus dem Dictionary
    rosen_base = TEST_FUNCTIONS["rosenbrock"]
    tf_rosen = fixed(rosen_base; n=2) # Wir fixieren sie auf 2D
    
    @testset "Rosenbrock 2D" begin
        res = analyze(tf_rosen)
        
        @test res isa SHGOResult
        # Diese Tests prüfen, ob das Objekt korrekt zurückgegeben wird
        @test res.num_basins >= 0
        @test length(res.local_minima) >= 0
        
        # Vergleich mit den Metadaten deiner Library
        @test min_value(tf_rosen) == 0.0
    end

    # 2. Test: Himmelblau (Festgelegte Dimension)
    himmel_base = TEST_FUNCTIONS["himmelblau"]
    tf_himmel = fixed(himmel_base) # Himmelblau ist fixiert (2D)
    
    @testset "Himmelblau" begin
        res = analyze(tf_himmel)
        
        @test res isa SHGOResult
        # Himmelblau hat laut Literatur 4 Minima
        @test name(tf_himmel) == "himmelblau"
    end
end
# End: test\test_pipeline.jl
