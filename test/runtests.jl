# File: test/runtests.jl
using SHGO
using Test

@testset "SHGO.jl Full Suite" begin
    # 1. Pipeline & Integration
    include("test_pipeline.jl")
       include("test_kuhn.jl") 
    # 2. Mathematischer Kern (Hardened)
    include("test_kuhn_hardened.jl")
    
   
  
	 include("test_heaps_validation.jl")
    include("shgogiar-test.jl")
    
end