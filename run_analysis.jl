# run_analysis.jl
# Starte mit: julia --project=. run_analysis.jl

using Pkg
Pkg.activate(".")

using SHGO
using NonlinearOptimizationTestFunctions
using Printf

println("=== SHGO Analyse - Six-Hump-Camelback (n=2) ===\n")

# Funktion laden
tf = fixed(TEST_FUNCTIONS["sixhumpcamelback"]; n=2)

# Parameter – du kannst hier rumspielen
n_div              = 15       # 8–15 ist ein guter Bereich (mehr → mehr Simplizes, langsamer)
use_grad_pruning   = false    # false = mehr Kandidaten, true = sehr aggressiv
verbose            = true
show_minima_detail = true     # Zeigt jedes gefundene lokale Minimum

println("Parameter:")
@printf "  Funktion:            %s\n" name(tf)
@printf "  Divisionen pro Achse: %d\n" n_div
@printf "  Gradient-Pruning:    %s\n" (use_grad_pruning ? "aktiviert" : "deaktiviert")
println("")

# Analyse starten
res = analyze(tf; 
    n_div = n_div,
    verbose = verbose,
    use_gradient_pruning = use_grad_pruning
)

# Ergebnis ausgeben
println("\n=== Ergebnis-Zusammenfassung ===")
@printf "Gefundene Basins:      %d\n" res.num_basins
@printf "Anzahl lokaler Minima: %d\n" length(res.local_minima)

if !isnothing(res.global_minimum)
    @printf "\nGlobales Minimum gefunden: f = %.6f\n" res.global_minimum.objective
    @printf "  Position: %.6f, %.6f\n" res.global_minimum.u[1] res.global_minimum.u[2]
else
    println("\nKein globales Minimum gefunden (noch zu wenige Kandidaten?)")
end

if show_minima_detail && !isempty(res.local_minima)
    println("\nGefundene lokale Minima (sortiert nach Objective):")
    sorted_minima = sort(res.local_minima, by = m -> m.objective)
    for (i, m) in enumerate(sorted_minima)
        @printf "  %2d.  f = %10.6f   @  (%.6f, %.6f)\n" i m.objective m.u[1] m.u[2]
    end
else
    println("\nKeine lokalen Minima gefunden – evtl. mehr Divisionen oder weniger Pruning probieren?")
end

println("\nFertig. Viel Erfolg beim Feintuning! 🦅")