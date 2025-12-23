# File: test/test_kuhn_hardened.jl (KORRIGIERTE VERSION)

@testset "Permutation Completeness & Uniqueness" begin
    # Test N=0 (edge case)
    @test isnothing(iterate(SHGO.KuhnPermutationIterator(0)))
    
    # Test N=1 through N=5
    for N in 1:5
        iter = SHGO.KuhnPermutationIterator(N)
        perms = collect(iter)
        
        # Korrekte Julia Syntax: Keine Strings nach der Bedingung!
        @test length(perms) == factorial(N) 
        
        # Alle müssen eindeutig sein
        @test length(unique(perms)) == factorial(N)
        
        # Validitäts-Check
        for perm in perms
            @test sort(perm) == SVector{N}(1:N)
        end
    end
end

@testset "Permutation Parity Balance" begin
    function count_inversions(perm)
        n = length(perm)
        inversions = 0
        for i in 1:n-1
            for j in i+1:n
                if perm[i] > perm[j]
                    inversions += 1
                end
            end
        end
        return inversions
    end
    
    for N in 2:4
        iter = SHGO.KuhnPermutationIterator(N)
        perms = collect(iter)
        
        parities = [count_inversions(p) % 2 for p in perms]
        n_even = count(==(0), parities)
        n_odd = count(==(1), parities)
        
        @test n_even == factorial(N) ÷ 2
        @test n_odd == factorial(N) ÷ 2
    end
end