# File: test\test_heaps_validation.jl
#test/test_heaps_validation.jl
using Test
using SHGO
using StaticArrays

function all_permutations_reference(n)
    if n == 1 return [Int[1]] end
    result = Vector{Vector{Int}}()
    for perm in all_permutations_reference(n-1)
        for i in 0:n-1
            new_p = copy(perm)
            insert!(new_p, i+1, n)
            push!(result, new_p)
        end
    end
    return result
end

@testset "Heap's Algorithm - Reference Validation" begin
    @testset "N=4: All Valid Permutations" begin
        expected = Set([SVector{4, Int}(p) for p in all_permutations_reference(4)])
        actual = Set(collect(SHGO.KuhnPermutationIterator(4)))
        @test actual == expected
        @test length(actual) == 24
    end

    @testset "Determinism & Edge Cases" begin
        iter = SHGO.KuhnPermutationIterator(3)
        @test iterate(iter)[1] == SVector(1, 2, 3)
        
        perms1 = collect(SHGO.KuhnPermutationIterator(4))
        perms2 = collect(SHGO.KuhnPermutationIterator(4))
        @test perms1 == perms2
    end
end
# End: test\test_heaps_validation.jl
