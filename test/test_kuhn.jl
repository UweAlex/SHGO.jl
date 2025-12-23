using Test
using SHGO
using StaticArrays
using NonlinearOptimizationTestFunctions

@testset "Kuhn Triangulation & Pruning" begin

    @testset "1. Geometrie: Permutationen (Heap's Algorithm)" begin
        # N=3 Check (sollten genau 6 sein)
        iter3 = SHGO.KuhnPermutationIterator(3)
        perms3 = collect(iter3)
        @test length(perms3) == 6
        @test length(unique(perms3)) == 6
        @test (@allocated collect(SHGO.KuhnPermutationIterator(3))) < 1000
    end

    @testset "2. Geometrie: Indizes" begin
        origin = CartesianIndex(1, 1)
        perm = SVector(2, 1) 
        indices = SHGO.generate_kuhn_indices(origin, perm)
        # Pfad: (1,1) -> (1,2) -> (2,2)
        @test indices == [CartesianIndex(1,1), CartesianIndex(1,2), CartesianIndex(2,2)]
    end

    @testset "3. Logik: Gradient Hull Pruning (Der 'scharfe' Test)" begin
        # Wir nutzen die Sphere-Funktion: f(x) = x^2 + y^2. Globales Min bei (0,0).
        # Gradient g(x) = 2x.
        
        # Setup: Sphere Funktion in 2D
        sphere_base = TEST_FUNCTIONS["sphere"]
        tf = fixed(sphere_base; n=2)
        
        # Fall A: Ein Bereich WEIT WEG vom Minimum (nur positive Gradienten)
        # Wir definieren Bounds so, dass wir weit im Positiven sind: [10, 12]
        # Ein Simplex hier hat Gradienten ~ [20, 20]. Die 0 ist NICHT in der Hülle.
        # -> Cache simulieren, der Indizes auf Koordinaten mappt
        # Wir hacken den Cache hier nicht, sondern nutzen SHGO's echten Cache.
        
        # Wir tricksen etwas mit den Bounds des Cache, um "weit weg" zu simulieren,
        # indem wir das Grid so definieren, dass Index (1,1) bei (10.0, 10.0) liegt.
        
        # Manuelles Cache-Setup für "Far Away" Szenario
        # LB = [10, 10], UB = [11, 11], Divisions = (1, 1)
        # Index (1,1) ist bei 10.0. Index (2,2) ist bei 11.0.
        cache_far = SHGO.VertexCache(
            Dict{CartesianIndex{2}, Tuple{Float64, SVector{2, Float64}}}(),
            ReentrantLock(),
            tf, 
            SVector(10.0, 10.0), # LB
            SVector(11.0, 11.0), # UB
            SVector(1.0, 1.0)    # cell_width
        )
        
        iter_pruned = SHGO.LazyKuhnSimplexes(CartesianIndex(1,1), cache_far)
        collected_far = collect(iter_pruned)
        
        # ERWARTUNG: Leere Liste! 
        # Weil alle Gradienten positiv sind, ist 0 nicht in der Hülle. Alles pruned.
        @test length(collected_far) == 0 
        
        # Fall B: Ein Bereich ÜBER dem Minimum
        # LB = [-0.5, -0.5], UB = [0.5, 0.5]
        # Das Gitter umschließt die 0. Gradienten zeigen in alle Richtungen.
        cache_hit = SHGO.VertexCache(
            Dict{CartesianIndex{2}, Tuple{Float64, SVector{2, Float64}}}(),
            ReentrantLock(),
            tf, 
            SVector(-0.5, -0.5), # LB
            SVector(0.5, 0.5),   # UB
            SVector(1.0, 1.0)    # cell_width
        )
        
        iter_keep = SHGO.LazyKuhnSimplexes(CartesianIndex(1,1), cache_hit)
        collected_hit = collect(iter_keep)
        
        # ERWARTUNG: Volle Anzahl Simplizes (2! = 2 Stück im 2D Rechteck)
        # Da 0 im Gradienten-Hull enthalten ist.
        @test length(collected_hit) == 2
        @test collected_hit[1] isa SHGO.Simplex
    end
end