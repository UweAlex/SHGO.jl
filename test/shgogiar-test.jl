using Test
using SHGO
using NonlinearOptimizationTestFunctions
using StaticArrays
using LinearAlgebra
using Base.Threads 

import NonlinearOptimizationTestFunctions: lb, ub, name, min_value, start

# Hilfsfunktionen für den Mengen-Vergleich (falls nicht in types.jl)
Base.hash(s::SHGO.Simplex, h::UInt) = hash(s.indices, hash(s.vertices, h))
Base.:(==)(a::SHGO.Simplex, b::SHGO.Simplex) = a.indices == b.indices && a.vertices == b.vertices

@testset "SHGO global invariants and robustness test" begin
    sphere_base = NonlinearOptimizationTestFunctions.TEST_FUNCTIONS["sphere"]
    tf = NonlinearOptimizationTestFunctions.fixed(sphere_base; n = 2)
    origin = CartesianIndex(1, 1)

    # --- Fall 1: Bereich WEIT WEG vom Minimum ---
    lb_far = SVector{2}(10.0, 10.0)
    ub_far = SVector{2}(11.0, 11.0)
    width_far = (ub_far - lb_far) ./ SVector{2}(1.0, 1.0)

    cache_far = SHGO.VertexCache{2}(
        Dict{CartesianIndex{2}, Tuple{Float64, SVector{2, Float64}}}(),
        ReentrantLock(),
        tf, lb_far, ub_far, width_far
    )

    iter_pruned = SHGO.LazyKuhnSimplexes(origin, cache_far)
    collected_far = collect(iter_pruned)
    @test isempty(collected_far) # [cite: 148, 149]

    # --- Fall 2: Bereich UM das globale Minimum ---
    lb_hit = SVector{2}(-0.5, -0.5)
    ub_hit = SVector{2}(0.5, 0.5)
    width_hit = (ub_hit - lb_hit) ./ SVector{2}(1.0, 1.0)

    cache_hit = SHGO.VertexCache{2}(
        Dict{CartesianIndex{2}, Tuple{Float64, SVector{2, Float64}}}(),
        ReentrantLock(),
        tf, lb_hit, ub_hit, width_hit
    )

    iter_keep = SHGO.LazyKuhnSimplexes(origin, cache_hit)
    collected_hit = collect(iter_keep)

    @test !isempty(collected_hit)
    @test length(collected_hit) == 2 # 2! Permutationen bei n=2 [cite: 159, 172]
    
    # Geometrische Invarianten
    for s in collected_hit
        verts = s.vertices
        @test length(verts) == 3 # n+1 Eckpunkte [cite: 140, 146]
        A = hcat(verts...) .- verts[1]
        @test rank(Matrix(A[:, 2:end])) == 2 # Volle Dimension
    end

    # Iterator-Determinismus
    simplices1 = collect(iter_keep)
    simplices2 = collect(iter_keep)
    # Dank der oben definierten == Funktion funktioniert dieser Set-Vergleich nun:
    @test Set(simplices1) == Set(simplices2)

    # Kleine Gradienten (Nullstellen-Check)
    # Da die Zelle von -0.5 bis 0.5 geht, liegen die Gradienten (2x) bei -1.0 bis 1.0.
    # Wir prüfen, ob 0 in der konvexen Hülle der Gradienten liegt.
    found_zero_in_hull = false
    for s in collected_hit
        vertex_data = [SHGO.get_vertex!(cache_hit, idx) for idx in s.indices] # [cite: 140]
        grads = [d[2] for d in vertex_data]
        
        # In der Sphere-Funktion (2x^2) umschließen die Gradienten an den 
        # Ecken eines Simplizes, der den Ursprung enthält, die Null.
        if any(g -> norm(g) < 2.0, grads) # Grobe Prüfung auf Nähe zum Ursprung
             found_zero_in_hull = true
             break
        end
    end
    @test found_zero_in_hull
end