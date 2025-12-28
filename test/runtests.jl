# File: test/runtests.jl
using Test  # <--- Das hat gefehlt!

@testset "Triangulation Full Suite" begin
    include("test_topical.jl")
    include("test_grid.jl")
end