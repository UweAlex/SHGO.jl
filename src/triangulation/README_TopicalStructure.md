
# README_TopicalStructure.md

# Dokumentation: TopicalStructure.jl

## 1. Übersicht

Die `TopicalStructure` ist das mathematische und datentechnische Fundament für die SHGO-Optimierung. Sie verwaltet die Diskretisierung des Suchraums als **Simplizialkomplex**. Sie fungiert als „Single Source of Truth“ (SSoT) zwischen der Geometrie (Raum), der Analysis (Energie) und der Topologie (Vernetzung). Im Gegensatz zu reinen Optimierern ist dieses System darauf ausgelegt, die gesamte Landschaft strukturell abzubilden.

## 2. Hierarchie der Objekte

### A. Vertex (Die Ecke)

Die kleinste Informationseinheit – das „Atom“ der Landschaft.

* **Eigenschaften:**
* `id`: Eindeutiger Primärschlüssel (entspricht der Grid-ID).
* `pos`: Statisch alloziertes Koordinaten-SVector (-D).
* `val`: Skalarer energetischer Zustand .


* **Rolle:**
* Startpunkt für lokale Optimierer.
* Referenzpunkt für die topologische Inzidenz (Star-Mapping).



### B. Simplex (Das Verknüpfungselement)

Eine topologische Einheit aus  Ecken – das „Molekül“.

* **Eigenschaften:**
* `id`: Laufende interne Kennung.
* `vertices`: Sortierte Liste der beteiligten Vertex-IDs.
* `min_val`: Cache des tiefsten Punktes im Simplex zur schnellen Filterung.


* **Rolle:**
* Definiert den lokalen Suchraum für die topologische Analyse.
* Dient als Container für die Identifikation vielversprechender Regionen (Kandidaten-Simplizes).



### C. TopicalManager (Das System)

Der zentrale Orchestrator, der alle Datenströme verwaltet.

* **Container:**
* `vertices`: `Dict{Int, Vertex}` – Das persistente Gedächtnis des Raums.
* `simplices`: `Dict{Int, Simplex}` – Die Landkarte der Vernetzung.
* `star_mapping`: `Dict{Int, Vector{Int}}` – Der **inverse Index** (Vertex  Simplizes). Dies ist die entscheidende Struktur für die Abfrage der Umgebung.
* `work_heap`: `PriorityQueue{Int, Float64}` – Die nach Tiefe sortierte To-Do-Liste der Vertices.


* **Fähigkeiten:**
* **Integritäts-Management**: Garantiert durch strikte Sortierung der IDs, dass Topologien eindeutig bleiben.
* **Nachbarschafts-Service**: Liefert über `get_star` sofort alle Flächen, die an einem Punkt hängen.



---

## 3. Dynamik & Phasen (Axiome)

Um die Korrektheit des SHGO-Algorithmus zu gewährleisten, folgt der Manager einem strikten Lebenszyklus. Die Trennung von Topologie-Aufbau und Basin-Identifikation ist hierbei zentral:

### Phase 1: Die Explorationsphase (Aufbau)

* **Aktion:** Gitterpunkte werden evaluiert und via `add_vertex!` registriert. Kuhn-Simplizes werden via `add_simplex!` gewebt.
* **Zustand:** Der `work_heap` füllt sich mit potenziellen Startpunkten. Die Vernetzung (Topologie) wird vollständig hergestellt.
* **Regel:** In dieser Phase findet noch keine lokale Optimierung statt. Wir bauen lediglich das "Skelett" der Landschaft.

### Phase 2: Die Analysephase (Konsum & Deduplizierung)

* **Aktion:** Der Algorithmus entnimmt die tiefste verfügbare Ecke aus dem `work_heap` (`consume_vertex!`).
* **Optimierung:** Von diesem Punkt aus startet eine lokale Suche.
* **Deduplizierung:** Das gefundene lokale Minimum wird mit bereits existierenden Minima verglichen. Führt die Suche zu einem neuen Punkt, wird ein neues **Attraction Basin** definiert.
* **Axiom:** Ein konsumierter Vertex dient als Keimzelle für eine lokale Suche. Mehrere Startpunkte können zum selben Basin führen – die Zusammenführung erfolgt über die Distanz im Phasenraum, nicht durch vorzeitiges topologisches Verschmelzen.

### Phase 3: Transparenz & Debugging

* **Offenheit:** Alle internen Datenstrukturen (Dicts) sind für Diagnosezwecke zugänglich, was eine Visualisierung der "Entdeckungsreise" ermöglicht.
* **Qualitäts-Sicherung:** Zukünftige Erweiterungen können ein `is_locked`-Flag nutzen, um Schreibzugriffe während der Analysephase zu verhindern.

---

## 4. Test-Szenarien (Qualm-Check)

Um die Stabilität der `TopicalStructure` zu garantieren, müssen folgende Tests bestanden werden:

1. **Inzidenz-Test:** Findet `get_star(v_id)` alle Simplizes, die diesen Vertex enthalten?
2. **Sortierungs-Test:** Erkennt das System, dass `[1, 2, 3]` und `[3, 2, 1]` denselben topologischen Körper beschreiben? (Vermeidung von Duplikaten).
3. **Heap-Integrität:** Liefert die `PriorityQueue` nach `consume_vertex!` wirklich den global tiefsten noch nicht bearbeiteten Punkt?

---

## 5. Vorbereitung zur Parallelisierung (Roadmap)

Obwohl die aktuelle Implementierung seriell arbeitet, ist die Architektur darauf ausgelegt, massiv parallele Workloads (z.B. auf 64 Kernen) ohne strukturellen Umbau zu unterstützen:

### A. Thread-Sicherheits-Strategie (Locking)

Da Julia-Dictionaries (`Dict`) nicht von Haus aus thread-sicher sind, bereitet das Design die Integration von Locks vor:

* **Granulares Locking:** Anstatt den gesamten Manager zu sperren, werden `ReentrantLocks` für das `star_mapping` und die Simplex-Tabellen vorgesehen.
* **Atomic Counters:** IDs (wie `_next_simplex_id`) werden als `Threads.Atomic{Int}` implementiert, um Race-Conditions bei der ID-Vergabe zu verhindern.

### B. Immutability (Unveränderlichkeit) der Atome

Die `Vertex`-Struktur ist als **immutable struct** definiert.

* **Vorteil:** Threads können Vertices lesen, ohne Memory-Barrieren oder Locks fürchten zu müssen, da die Daten (`pos`, `val`) nach der Erstellung garantiert unveränderlich sind. Dies ist die Basis für "Lock-free Reading".

### C. Partitionierung des Work-Heaps

Der globale `work_heap` kann für parallele Suchen in regionale oder Thread-lokale Heaps partitioniert werden. Die `consume_vertex!`-Logik ist so gekapselt, dass sie später durch eine Lastverteilungs-Strategie ersetzt werden kann.

### D. Batch-Insertion (Kuhn-Parallelisierung)

Kuhn-Simplizes können pro Gitterzelle völlig unabhängig berechnet werden.

* **Vorbereitung:** Threads sammeln ihre generierten Simplizes lokal und führen am Ende einer Epoche eine "Batch-Insertion" in den `TopicalManager` durch. Dies minimiert die Zeit, in der globale Locks gehalten werden müssen.

---

**Abschluss-Statement zur Qualifizierung:**
Durch die Verwendung von `StaticArrays` (SVector) und dem strikten Verzicht auf globale Zustände innerhalb der Kernfunktionen ist der Code bereits "pure" genug, um unmittelbar in `Threads.@threads`-Schleifen genutzt zu werden, sobald die Schreibzugriffe auf die zentralen Register abgesichert sind.