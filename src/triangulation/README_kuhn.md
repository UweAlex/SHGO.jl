
# README: Kuhn-Triangulierung & Simplex-Generierung

## 1. Status Quo (`kuhn.jl`)

Die Implementierung ist eine leistungsfähige, in sich geschlossene Einheit, die eine Hyperkubus-Zelle in  Simplizes zerlegt. Sie ist vollständig in das neue Framework integriert und operiert auf der `GridStructure`.

### Kern-Komponenten:

* **`KuhnPermutationIterator{N}`**: Erzeugt statisch alle Permutationen der Einheitsvektoren, um Allokationen während der Laufzeit zu minimieren. Dies ist entscheidend für die Performance in höheren Dimensionen.
* **`generate_kuhn_simplices`**: Die Hauptfunktion, die das Gitter durchläuft und für jede Zelle die topologischen Pfade berechnet. Sie ist nun entkoppelt von der Werterhebung und arbeitet rein auf Basis von IDs und Geometrie-Vorgaben.

---

## 2. Das "Gradienten-Abrakadabra" (Pruning)

Ein besonderes Merkmal von SHGO.jl ist das integrierte **Gradient Hull Pruning**. Dies ist eine Erweiterung gegenüber dem Standard-SHGO, die die Effizienz massiv steigern kann.

### Das Prinzip:

Mathematisch gesehen muss an einem lokalen Minimum der Gradient Null sein. Das Pruning nutzt die konvexe Hülle der Gradienten an den Ecken eines Simplex:

1. Es werden die Gradienten  der Simplex-Ecken erhoben.
2. Es wird geprüft, ob der Ursprung (Nullpunkt) innerhalb dieser konvexen Hülle liegt.
3. **Die Logik:** Wenn alle Gradienten in eine ähnliche Richtung zeigen (die Hülle also die Null nicht enthält), "rutscht" die Funktion in diesem Simplex nur in eine Richtung ab. Ein lokales Minimum ist dort topologisch unmöglich. Der Simplex wird verworfen.

### Problem & Lösung:

* **Problem:** Das Pruning kann bei komplexen Landschaften (wie dem Six-Hump-Camelback) zu aggressiv sein, wenn die Auflösung (`n_div`) zu gering ist. Kleine, schmale Basins können "verschluckt" werden, wenn die Gradienten an den Ecken des Simplex den Nullpunkt knapp verfehlen.
* **Status:** Das Pruning ist über den Parameter `use_gradient_pruning::Bool` optional steuerbar. In der Diagnose-Phase (Six-Hump) wird es oft deaktiviert, um die vollständige strukturelle Abdeckung zu garantieren.

---

## 3. Die Architektur von `Kuhn.jl`

Kuhn fungiert im Gesamtsystem als **reiner topologischer Dienstleister** (Worker). Die Architektur folgt diesen Prinzipien:

1. **ID-basierte Operation:** Kuhn berechnet keine physikalischen Positionen mehr. Er nutzt das `GridStructure`-Objekt, um mittels `get_linear_id` die global eindeutigen Identifikatoren der Vertices zu bestimmen.
2. **Direkte Meldung:** Kuhn erzeugt keine großen temporären Arrays von Simplex-Objekten mehr. Jeder valide gefundene Simplex wird direkt an den `TopicalManager` via `add_simplex!(tm, ids)` gemeldet.
3. **Zustandslosigkeit:** Das Modul speichert keine Daten über den Zustand der Welt; es "webt" lediglich die Verbindungen zwischen den bereits existierenden Vertices.

---

## 4. Implementierter Flow

Die Zusammenarbeit zwischen den Modulen ist nun wie folgt festgeschrieben:

* **Schritt 1:** Der Orchestrator (`analyze`) lässt das `Grid` die Welt besiedeln (`add_vertex!`).
* **Schritt 2:** Kuhn wird pro Gitterzelle aufgerufen.
* **Schritt 3:** Kuhn identifiziert die  Simplizes einer Zelle.
* **Schritt 4 (Optional):** Falls `use_gradient_pruning` aktiv ist, werden nur Simplizes an den `TopicalManager` übertragen, deren Gradienten-Hülle den Nullpunkt enthält.
* **Schritt 5:** Der `TopicalManager` speichert die validen Simplizes und aktualisiert den `star_mapping` (inverser Index).

---

### Qualifizierte Bewertung:

Die Migration vom "Legacy-Kern" zum neuen Framework ist abgeschlossen. Die Entkoppelung ermöglicht es nun, verschiedene Triangulierungs-Verfahren (nicht nur Kuhn) einzusetzen, ohne den Rest des Systems (den Manager oder den Solver) anpassen zu müssen. Das Gradient-Pruning bleibt als mächtiges Werkzeug erhalten, erfordert jedoch eine sorgfältige Abstimmung mit der Gitter-Auflösung.

