# File: test/test_grid.jl
using Test
using StaticArrays

const GRID_FILE = joinpath(@__DIR__, "..", "src", "triangulation", "Grid.jl")
include(GRID_FILE)
using .Grid

@testset "Grid.jl - Brutal Testing" begin

    @testset "1. Numerische Präzision (Extreme Skalen)" begin
        # Teste winzige Abstände (Nanometer-Bereich)
        lower = [0.0, 0.0]
        upper = [1e-9, 1e-9]
        # Um eine Schrittweite von exakt 1e-10 zu erhalten, 
        # brauchen wir bei einer Spanne von 1e-9 genau 10 Intervalle, also 11 Punkte.
        n_points = 11 
        dims = [n_points, n_points]
        grid = GridStructure(lower, upper, dims)
        
        pos_end = get_vertex_pos(grid, [n_points, n_points])
        @test all(pos_end .≈ 1e-9)
        @test grid.steps[1] ≈ 1e-10
    end

    @testset "2. Out-of-Bounds Angriffe" begin
        grid = GridStructure([0.0], [1.0], [5])
        
        # Teste Index 0 (Julia ist 1-basiert)
        @test_throws BoundsError get_vertex_pos(grid, [0])
        # Teste Index n+1
        @test_throws BoundsError get_vertex_pos(grid, [6])
        # Teste negative Indizes
        @test_throws BoundsError get_vertex_pos(grid, [-1])
    end

    @testset "3. Degenerierte Gitter" begin
        # Was passiert, wenn upper < lower?
        @test_throws ArgumentError GridStructure([1.0], [0.0], [5])
        
        # Was passiert bei nur einem Punkt pro Dimension? (Jetzt durch Check im Konstruktor abgefangen)
        @test_throws ArgumentError GridStructure([0.0], [1.0], [1])
    end

    @testset "4. Massive Dimensionen (Kombinatorik)" begin
        # 10D Gitter mit nur 2 Punkten pro Achse
        N = 10
        grid = GridStructure(zeros(N), ones(N), fill(2, N))
        @test get_total_vertices(grid) == 2^10 # 1024
        
        # Checke die "letzte Ecke" in 10D
        last_corner = get_vertex_pos(grid, fill(2, N))
        @test all(last_corner .== 1.0)
    end
end