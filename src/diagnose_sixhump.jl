# File: diagnose_sixhump.jl
using Pkg
Pkg.activate(".")
using SHGO, NonlinearOptimizationTestFunctions, Printf

tf = fixed(TEST_FUNCTIONS["sixhumpcamelback"]; n=2)

println("="^80)
println("SIX-HUMP-CAMELBACK DIAGNOSE")
println("="^80)

# Bekannte Minima aus der Literatur zur Validierung
known_minima = [([-0.0898, 0.7126], -1.0316), ([0.0898, -0.7126], -1.0316)]

for n_div in [10, 15]
    println("\n--- n_div = $n_div ---")
    res = analyze(tf; n_div = n_div)
    @printf "Gefundene Basins: %d\n" res.num_basins
end
# End: diagnose_sixhump.jl