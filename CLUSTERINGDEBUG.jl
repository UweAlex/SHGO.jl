using Pkg
Pkg.activate(".")

using SHGO
using NonlinearOptimizationTestFunctions
using Printf

println("="^80)
println("CLUSTERING DEBUG - Six-Hump-Camelback")
println("="^80)

tf = fixed(TEST_FUNCTIONS["sixhumpcamelback"]; n=2)

# Test mit moderater Auflösung
n_div = 12
println("\n--- Analyse mit n_div = $n_div, OHNE Gradient-Pruning ---\n")

# Führe Analyse durch
res = analyze(tf; n_div = n_div, verbose = true, use_gradient_pruning = false)

println("\n" * "="^80)
println("ERGEBNISSE")
println("="^80)
@printf "Gefundene Basins:    %d\n" res.num_basins
@printf "Lokale Minima:       %d\n" length(res.local_minima)

if !isempty(res.local_minima)
    println("\nGefundene Minima (nach Objective sortiert):")
    sorted = sort(res.local_minima, by = m -> m.objective)
    for (i, m) in enumerate(sorted)
        @printf "  %2d. f = %8.4f  @ (%.4f, %.4f)\n" i m.objective m.u[1] m.u[2]
    end
    
    # Prüfe Unique-Werte (vielleicht Duplikate?)
    unique_objectives = unique([m.objective for m in res.local_minima])
    @printf "\nAnzahl UNIQUE Objective-Werte: %d\n" length(unique_objectives)
    
    if length(unique_objectives) < length(res.local_minima)
        println("⚠️  WARNUNG: Es gibt Duplikate! Clustering findet mehrfach dasselbe Minimum.")
    end
end

println("\n" * "="^80)
println("DIAGNOSE")
println("="^80)
println("Erwartung:  6 verschiedene Minima")
println("Realität:   Siehe oben")
println("\nMögliche Probleme:")
println("  1. Gradient-Pruning wirft zu viele Simplizes weg")
println("  2. Clustering verbindet alles zu einem Basin")
println("  3. Grid zu grob → lokale Minima nicht erfasst")
println("  4. Lokale Optimierung konvergiert mehrfach zum selben Punkt")
println("="^80)