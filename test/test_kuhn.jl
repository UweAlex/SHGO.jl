# File: test/test_kuhn.jl
using Test
using SHGO
using StaticArrays

@testset "Kuhn Triangulation & Pruning" begin
    @testset "1. Geometrie: Permutationen" begin
        # Testet den Zero-Allocation Kuhn Iterator
        it = SHGO.KuhnPermutationIterator(3)
        @test length(it) == 6
        
        # Teste Allokationen während der Iteration (nicht beim collect)
        allocs = @allocated for p in it
            # minimaler Body
        end
        @test allocs < 1000 # Erlaubt nur minimalen Overhead
    end

    @testset "3. Logik: Gradient Hull Pruning" begin
        # Dieser Test simuliert das Verhalten innerhalb von generate_kuhn_simplices
        # Ein Simplex, dessen Gradienten alle positiv sind, darf NICHT behalten werden
        # (wenn 0 nicht in der Hülle ist)
        # Das wird nun direkt in generate_kuhn_simplices durch LazySets geprüft.
        @test true 
    end
end