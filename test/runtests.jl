#test/runtests.jl

using SHGO
using Test

@testset "SHGO.jl" begin
    include("test_pipeline.jl")
end