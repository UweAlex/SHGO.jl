# File: maintain.jl
using Dates

function ultra_safe_maintenance(root_dir = pwd())
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    log_file = joinpath(root_dir, "maintenance_$timestamp.log")
    
    open(log_file, "w") do log
        write(log, "Maintenance Start: $(now())\n" * "-"^40 * "\n")
        
        for (root, dirs, files) in walkdir(root_dir)
            # Versteckte Ordner wie .git ignorieren
            if occursin(".git", root) continue end
            
            for file in files
                if endswith(file, ".jl") && file != "maintenance_script.jl"
                    full_path = joinpath(root, file)
                    rel_path = relpath(full_path, root_dir)
                    
                    try
                        process_file_ultra_safe(full_path, rel_path, timestamp, log)
                    catch e
                        write(log, "ERROR bei $rel_path: $e\n")
                    end
                end
            end
        end
    end
    println("Abgeschlossen. Details findest du in: $log_file")
end

function process_file_ultra_safe(full_path, rel_path, ts, log_io)
    lines = readlines(full_path, keep=true)
    orig_count = length(lines)
    modified = false
    
    header = "# File: $rel_path\n"
    footer = "# End: $rel_path\n"

    # Header-Check (erste 5 Zeilen)
    if !any(occursin("# File: $rel_path", l) for l in lines[1:min(5, end)])
        insert!(lines, 1, header)
        modified = true
    end

    # Footer-Check (letzte 5 Zeilen)
    if !any(occursin("# End: $rel_path", l) for l in lines[max(1, end-4):end])
        if !isempty(lines) && !endswith(lines[end], "\n")
            lines[end] *= "\n"
        end
        push!(lines, footer)
        modified = true
    end

    if modified
        backup_path = full_path * ".$ts.bak"
        
        # 1. Backup schreiben
        cp(full_path, backup_path, force=false)
        
        if isfile(backup_path)
            # 2. Original schreiben
            open(full_path, "w") do f
                foreach(l -> write(f, l), lines)
            end
            
            # 3. Verifizieren
            new_count = length(readlines(full_path))
            if new_count >= orig_count
                write(log_io, "SUCCESS: $rel_path (Backup: $(basename(backup_path)))\n")
                println("✓ $rel_path")
            else
                write(log_io, "WARNING: $rel_path könnte korrupt sein (Zeilen geschrumpft!)\n")
            end
        end
    else
        write(log_io, "SKIP: $rel_path (Bereits aktuell)\n")
    end
end

ultra_safe_maintenance()
# End: maintain.jl
