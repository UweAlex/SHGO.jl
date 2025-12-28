
# SHGO.jl

**Stochastic Homology Global Optimization in Julia**

SHGO.jl ist nicht nur ein globaler Optimierer, sondern ein Instrument zur **Strukturanalyse von Optimierungslandschaften**. Es basiert auf der Idee, die topologische Information einer Zielfunktion zu nutzen, um alle lokalen Minima effizient zu finden und das globale Minimum zu garantieren.

---

## ⚠️ Der "Heureka-Moment": Korrektur der Basin-Theorie

In der frühen Entwicklung gab es einen fundamentalen logischen Fehler beim Clustering. Diese Erkenntnis ist der Kern der aktuellen Architektur:

### Das Problem: Premature Clustering

Früher wurde versucht, **aktive Simplizes** (Simplizes, in denen ein Minimum vermutet wird) zu gruppieren, *bevor* die lokale Optimierung stattfand.

* **Fehler:** "Multiplikative Barrieren". Zwei Simplizes können topologisch benachbart sein, aber zu zwei völlig verschiedenen Minima führen.
* **Folge:** Basins wurden "verschmolzen", bevor sie existierten. Die Anzahl der gefundenen Minima war instabil.

### Die Lösung: Der SHGO-Flow (Korrekt)

1. **Sampling:** Den Raum in ein Gitter aus Simplizes unterteilen.
2. **Topologische Filterung:** Identifikation "aktiver" Simplizes (Kandidatenregionen).
3. **Lokale Optimierung:** *Jeder* aktive Simplex startet eine lokale Suche.
4. **Deduplizierung (Echtes Clustering):** Erst die *Resultate* der Optimierung werden anhand ihrer Distanz im Phasenraum gruppiert.

---

## 📘 Terminologie & Taxonomie

Um Missverständnisse zu vermeiden (insbesondere im Vergleich zur SciPy-Implementierung), nutzt SHGO.jl folgende Definitionen:

| Begriff | Definition | Rolle im Algorithmus |
| --- | --- | --- |
| **Simplex** | Kleinste geometrische Einheit des Gitters. | Datenspeicher (Werte/Gradienten). |
| **Star-Domain** | Die Nachbarschaft um einen Vertex. | Basis für die Homologie-Analyse. |
| **Kandidaten-Region** | Zusammenhängende Menge aktiver Simplizes. | Ein "topologisches Feature" der Landschaft. |
| **Attraction Basin** | Menge aller Punkte, die zum selben Minimum konvergieren. | **Wird erst nach Schritt 4 (Deduplizierung) gezählt.** |

---

## 🛠 Architektur & Design-Prinzipien

### 1. Separation of Concerns

* **TopicalManager:** Verwaltet die Geometrie und Topologie (Vertices, IDs, Simplex-Beziehungen).
* **Solver-Abstraktion:** Die lokale Optimierung ist entkoppelt.
* **Analyse-Layer:** Liefert ein "Profil" der Landschaft, nicht nur eine Zahl.

### 2. Eager vs. Lazy Evaluation

Im Gegensatz zu reinen Optimierern (die "Lazy" arbeiten, um Rechenzeit zu sparen), verfolgt SHGO.jl oft einen **Eager-Ansatz**:

* **Ziel:** Vollständige Abdeckung der Landschaft.
* **Vorteil:** Wir erhalten eine verlässliche Verteilung der Funktionswerte, was für die Diagnose von "Deceptiveness" (Irreführung) der Funktion kritisch ist.

### 3. Julia-spezifische Optimierungen

* **StaticArrays:** Für blitzschnelle geometrische Berechnungen in niedrigen Dimensionen.
* **Type-Safety:** Klare Trennung zwischen Vertex-IDs und physischen Koordinaten zur Vermeidung von Floating-Point-Fehlern in der Topologie-Logik.

---

## 🚀 Status Quo (Six-Hump-Camelback Test)

Die aktuelle Version löst das klassische **Six-Hump-Camelback** Problem (2D) absolut stabil:

* **Erwartet:** 6 lokale Minima (davon 2 global).
* **Resultat:** 6/6 Basins werden verlässlich gefunden (`res.num_basins == 6`).
* **Differenzierung:** Der Algorithmus erkennt die Struktur auch bei geringer Auflösung (`n_div = 12`), ohne Minima zu übersehen oder künstlich aufzublähen.

---

## 📈 Roadmap

* [x] Korrektes Basin-Clustering (Post-Optimization).
* [ ] Adaptive Gitter-Verfeinerung (Local Refinement).
* [ ] Integration von `ForwardDiff` für exakte lokale Gradienten.
* [ ] Parallelisierung der lokalen Suchen.

---

### Warum SHGO.jl?

*SciPy will lösen. SHGO.jl will verstehen.* Dieses Projekt ist für Anwender gedacht, die nicht nur wissen wollen, *wo* das Minimum liegt, sondern *wie* die gesamte energetische Landschaft beschaffen ist.