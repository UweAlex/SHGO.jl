module Basins
using Graphs

export find_basins

"""
    find_basins(tm, active_s_ids, N_dim)

Clustert Simplizes zu Basins mit ADAPTIVER Barrieren-Erkennung.

STRATEGIE:
1. Berechne globale Value-Range der Funktion
2. Setze Toleranz relativ zu dieser Range (z.B. 5%)
3. Zwei Simplizes gehören zum selben Basin, wenn:
   - Sie topologisch benachbart sind (N gemeinsame Vertices)
   - Die Barriere zwischen ihnen "klein" ist relativ zur globalen Landschaft
"""
function find_basins(tm, active_s_ids::Vector{Int}, N_dim::Int)
    n = length(active_s_ids)
    n == 0 && return Vector{Vector{Int}}()
    
    # 1. Berechne globale Value-Range für adaptive Toleranz
    all_values = Float64[]
    for s_id in active_s_ids
        s = tm.simplices[s_id]
        for vid in s.vertices
            push!(all_values, tm.vertices[vid].val)
        end
    end
    
    f_min_global = minimum(all_values)
    f_max_global = maximum(all_values)
    value_range = f_max_global - f_min_global
    
    # KRITISCHER PARAMETER: Adaptive Toleranz basierend auf der Landschaft
    # Für Six-Hump: Range ≈ 3.1 (von -1.03 bis +2.1)
    # → 5% davon = 0.155 (sollte globale von lokalen Minima trennen)
    relative_tolerance = 0.05  # 5% der Range
    barrier_tolerance = value_range * relative_tolerance
    
    g = Graphs.SimpleGraph(n)
    
    for i in 1:n, j in (i+1):n
        s1 = tm.simplices[active_s_ids[i]]
        s2 = tm.simplices[active_s_ids[j]]
        
        shared = intersect(s1.vertices, s2.vertices)
        
        # Topologischer Check: Müssen eine Facette teilen
        if length(shared) >= N_dim
            # Energetischer Check: Wie hoch ist die Barriere?
            f_interface = minimum(tm.vertices[v].val for v in shared)
            f_min1 = s1.min_val
            f_min2 = s2.min_val
            
            # Barrieren-Höhe: Wie viel höher ist die Grenze als das höhere Minimum?
            barrier_height = f_interface - max(f_min1, f_min2)
            
            # Wenn die Barriere klein genug ist → selbes Basin
            if barrier_height <= barrier_tolerance
                Graphs.add_edge!(g, i, j)
            end
        end
    end
    
    comps = Graphs.connected_components(g)
    return [[active_s_ids[idx] for idx in c] for c in comps]
end

end # module