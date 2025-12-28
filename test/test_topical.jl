using Test
using StaticArrays
using DataStructures

const TOPICAL_FILE = joinpath(@__DIR__, "..", "src", "triangulation", "TopicalStructure.jl")
include(TOPICAL_FILE)
using .TopicalStructure

@testset "TopicalStructure - Deep Dive Stress Test" begin

    @testset "1. ID-Eindeutigkeit & Stabilität" begin
        tm = TopicalManager{3}()
        ids = Set{Int}()
        for i in 1:100
            v = add_vertex!(tm, SVector(rand(3)...), rand())
            push!(ids, v.id)
        end
        @test length(ids) == 100  # Keine Dubletten
        @test maximum(ids) == 100 # Fortlaufende Zählung
    end

    @testset "2. Topologische Nachbarschaft (2D Netz)" begin
        # Wir bauen manuell ein kleines Netz aus 4 Simplizes (Quadrat aus 2 Dreiecken)
        tm = TopicalManager{2}()
        # Vertices für ein 2x2 Gitter
        v1 = add_vertex!(tm, SVector(0.0, 0.0), 0.1).id
        v2 = add_vertex!(tm, SVector(1.0, 0.0), 0.2).id
        v3 = add_vertex!(tm, SVector(0.0, 1.0), 0.3).id
        v4 = add_vertex!(tm, SVector(1.0, 1.0), 0.4).id
        
        s1 = add_simplex!(tm, [v1, v2, v3])
        s2 = add_simplex!(tm, [v2, v3, v4])
        
        # Ein Vertex in der Mitte (v2 oder v3) muss in beiden Simplizes sein
        star_v2 = get_star(tm, v2)
        @test length(star_v2) == 2
        @test s1 in star_v2 && s2 in star_v2
        
        # Ein Eck-Vertex (v1) darf nur in einem sein
        @test length(get_star(tm, v1)) == 1
    end

    @testset "3. Heap-Stress: Massives Consuming" begin
        tm = TopicalManager{1}()
        n = 500
        # Füge 500 Vertices mit zufälligen Werten hinzu
        for i in 1:n
            add_vertex!(tm, SVector(Float64(i)), rand() * 100.0)
        end
        
        last_val = -1.0
        for i in 1:n
            v_id = next_work_vertex(tm)
            current_val = tm.vertices[v_id].val
            @test current_val >= last_val # Die Werte MÜSSEN aufsteigend aus dem Heap kommen
            last_val = current_val
            consume_vertex!(tm, v_id)
        end
        @test next_work_vertex(tm) === nothing # Heap muss leer sein
    end

    @testset "4. Robustheit: Illegale Simplizes" begin
        tm = TopicalManager{3}() # 3D braucht 4 Vertices
        v1 = add_vertex!(tm, SVector(0.,0.,0.), 1.)
        v2 = add_vertex!(tm, SVector(1.,0.,0.), 1.)
        v3 = add_vertex!(tm, SVector(0.,1.,0.), 1.)
        
        # Teste Unter-Dimensionierung
        @test_throws ArgumentError add_simplex!(tm, [v1.id, v2.id, v3.id])
        
        # Teste doppelte Vertices im selben Simplex (entarteter Simplex)
        # Manche Algorithmen erlauben das, aber für SHGO ist es oft ein Fehler
        # Hier prüfen wir, ob dein Manager das (noch) zulässt oder ob wir eine 
        # Validierung einbauen wollen:
        @test_throws ArgumentError add_simplex!(tm, [v1.id, v1.id, v2.id, v3.id]) 
        # Hinweis: Falls dieser Test fehlschlägt, müssen wir die Logik in add_simplex! 
        # um 'length(unique(vertex_ids)) == N+1' erweitern.
    end

    @testset "5. Speicher-Integrität nach Löschung" begin
        tm = TopicalManager{2}()
        v1 = add_vertex!(tm, SVector(0.0, 0.0), 1.0)
        consume_vertex!(tm, v1.id)
        @test is_consumed(tm, v1.id)
        # Sicherstellen, dass der Vertex noch im Dict ist, aber nicht mehr im Heap
        @test haskey(tm.vertices, v1.id)
        @test isempty(tm.work_heap)
    end
end