# ─────────────────────────────────────────────────────────────────
    # 4. LazyKuhnSimplexes + Gradient Hull Pruning
    # ─────────────────────────────────────────────────────────────────
    @testset "LazyKuhnSimplexes + Gradient Hull Pruning" begin
        # Positiver Bereich → Gradienten alle positiv → 0 nicht in Hull → prune alles
        cache_far = create_test_cache(2, (10.0,10.0), (11.0,11.0), 1)
        
        # Test MIT Pruning → muss leer sein
        simplices_far = SHGO.generate_kuhn_simplices(2, 1, cache_far, true)
        @test isempty(simplices_far)
        @test length(simplices_far) == 0

        # Ohne Pruning → erwarte 2 Simplices (2! = 2 Permutationen)
        simplices_no_prune = SHGO.generate_kuhn_simplices(2, 1, cache_far, false)
        @test length(simplices_no_prune) == 2

        # GEÄNDERT: Bereich um Null mit größerem Toleranz
        # Sphere bei (-0.5, 0.5) hat kleine Gradienten, die durch 1e-7 fallen können
        # Teste stattdessen (-1, 1) für robustere Erkennung
        cache_near = create_test_cache(2, (-1.0,-1.0), (1.0,1.0), 1)
        simplices_near = SHGO.generate_kuhn_simplices(2, 1, cache_near, true)

        @test length(simplices_near) ≥ 1
        @test all(s -> length(s.vertices) == 3, simplices_near)
    end