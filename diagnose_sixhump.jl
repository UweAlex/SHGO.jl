# diagnose_sixhump.jl
# Detaillierte Diagnose der Six-Hump-Camelback Probleme

using Pkg
Pkg.activate(".")

using SHGO
using NonlinearOptimizationTestFunctions
using Printf
using Statistics

println("="^80)
println("SIX-HUMP-CAMELBACK DIAGNOSE")
println("="^80)

tf = fixed(TEST_FUNCTIONS["sixhumpcamelback"]; n=2)

# Bekannte Minima aus der Literatur
known_minima = [
    ([-0.0898, 0.7126], -1.0316),
    ([0.0898, -0.7126], -1.0316),
    ([-1.7036, 0.7961], -0.2155),
    ([1.7036, -0.7961], -0.2155),
    ([-1.6071, -0.5687], 2.1040),
    ([1.6071, 0.5687], 2.1040)
]

println("\nBekannte Minima:")
for (i, (pos, val)) in enumerate(known_minima)
    @printf "  %d. f(%.4f, %.4f) = %.4f\n" i pos[1] pos[2] val
end

println("\n" * "="^80)
println("TEST 1: Verschiedene Grid-Auflösungen (OHNE Refinement) - FIXED VERSION")
println("="^80)

for n_div in [5, 8, 10, 12, 15,20,30]
    println("\n--- n_div = $n_div ---")
    
    res = analyze(tf; 
        n_div = n_div,
        verbose = false,
        use_gradient_pruning = false,
        refinement_levels = 0  # Erstmal ohne
    )
    
    @printf "Simplizes insgesamt: ~%d (theoretisch: %d Zellen * 2! = %d)\n" (n_div^2 * 2) (n_div^2) (n_div^2 * 2)
    @printf "Gefundene Basins:    %d\n" res.num_basins
    @printf "Lokale Minima:       %d\n" length(res.local_minima)
    
    if !isempty(res.local_minima)
        println("\nGefundene Minima (sortiert):")
        sorted = sort(res.local_minima, by = m -> m.objective)
        for (i, m) in enumerate(sorted)
            @printf "  %d. f = %8.4f  @ (%.4f, %.4f)\n" i m.objective m.u[1] m.u[2]
        end
        
        # Pruefe Abstand zu bekannten Minima
        global_found = abs(sorted[1].objective + 1.0316) < 0.01
        println("\nGlobales Minimum gefunden: $(global_found ? "ja" : "nein")")
    end
end

println("\n" * "="^80)
println("TEST 2: Mit Gradient-Pruning")
println("="^80)

for use_pruning in [false, true]
    println("\n--- Gradient-Pruning: $use_pruning ---")
    
    res = analyze(tf; 
        n_div = 12,
        verbose = false,
        use_gradient_pruning = use_pruning,
        refinement_levels = 0
    )
    
    @printf "Gefundene Basins:  %d\n" res.num_basins
    @printf "Lokale Minima:     %d\n" length(res.local_minima)
    
    if !isempty(res.local_minima)
        best = minimum(m.objective for m in res.local_minima)
        @printf "Bestes Minimum:    %.6f\n" best
    end
end

println("\n" * "="^80)
println("TEST 3: Mit Refinement (kann langsam sein!)")
println("="^80)

for ref_level in [0, 1, 2]
    println("\n--- Refinement Level: $ref_level ---")
    
    res = analyze(tf; 
        n_div = 8,  # Kleiner, weil Refinement multipliziert
        verbose = false,
        use_gradient_pruning = false,
        refinement_levels = ref_level
    )
    
    @printf "Gefundene Basins:  %d\n" res.num_basins
    @printf "Lokale Minima:     %d\n" length(res.local_minima)
    
    if !isempty(res.local_minima)
        best = minimum(m.objective for m in res.local_minima)
        worst = maximum(m.objective for m in res.local_minima)
        @printf "Range:             %.4f bis %.4f\n" best worst
    end
end

println("\n" * "="^80)
println("TEST 4: Clustering-Detail-Analyse")
println("="^80)

# Manuell ein paar Simplizes erzeugen und Clustering testen
println("\nGeneriere Simplizes für n_div=10...")
res = analyze(tf; n_div=10, verbose=true, refinement_levels=0)

println("\n" * "="^80)
println("ZUSAMMENFASSUNG")
println("="^80)
println("Erwartung:  6 Basins (aus Literatur)")
println("Realität:   Siehe oben - vermutlich deutlich weniger")
println("\nMögliche Ursachen:")
println("  1. share_face() zu restriktiv/zu permissiv")
println("  2. Value-based Pruning zu aggressiv")
println("  3. Lokale Optimierung konvergiert nicht richtig")
println("  4. Grid zu grob für kleine Basins")
println("="^80)