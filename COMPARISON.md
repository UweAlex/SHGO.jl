

# Vergleich: SHGO.jl vs. SciPy SHGO (Python)

Dieses Dokument dient dazu, die architektonischen und konzeptionellen Unterschiede zwischen der ursprünglichen SciPy-Implementierung (`scipy.optimize.shgo`) und der Julia-Version `SHGO.jl` darzulegen.

## 1. Philosophischer Ansatz: Lösen vs. Verstehen

* **SciPy SHGO** ist als klassischer **Optimierer** konzipiert. Das Ziel ist es, so schnell wie möglich das globale Minimum zu finden. Sobald die homologische Abbruchbedingung erfüllt ist, stoppt der Prozess.
* **SHGO.jl** versteht sich als **Landschafts-Analysator**. Die Prämisse ist, dass die Struktur der Funktion (alle lokalen Minima, deren Einzugsgebiete und die topologische Konnektivität) für den Anwender oft wertvoller ist als nur der tiefste Punkt.

## 2. Technischer Vergleich

| Feature | SciPy SHGO (Python) | SHGO.jl (Julia) |
| --- | --- | --- |
| **Sprachbasis** | Python / C (Monolithisch) | Pure Julia (Modular) |
| **Datenstrukturen** | Schwerfällige Python-Objekte für Vertices und Simplizes. | **ID-basierte Topologie** mit **StaticArrays** für maximale Cache-Effizienz. |
| **Evaluation** | **Lazy Evaluation**: Berechnet nur, was für den nächsten Schritt nötig ist. | **Eager Evaluation**: Erstellt ein vollständiges Gitter-Profil zur Strukturdiagnose. |
| **Clustering** | **In-Graph Clustering**: Versucht Simplizes während des Wachstums zu gruppieren (Fehleranfällig bei komplexen Barrieren). | **Post-Optimization Clustering**: Nutzt die Resultate der lokalen Suche zur exakten Basin-Definition. |
| **Parallelisierung** | Durch den Global Interpreter Lock (GIL) limitiert; Parallelisierung oft auf Solver-Ebene "aufgepfropft". | Nativ auf **Multithreading** ausgelegt. Jeder aktive Simplex kann unabhängig optimiert werden. |
| **Erweiterbarkeit** | Schwer erweiterbar, da Topologie und Algorithmus eng verkoppelt sind. | Modularer `TopicalManager` erlaubt das Austauschen von Sampling-Strategien und Solvern. |

## 3. Die Basin-Definition (Der "Core-Fix")

Der wichtigste Fortschritt in `SHGO.jl` betrifft die Behandlung von **Attraction Basins**:

In der Python-Version führt eine hohe `n_div` (Auflösung) oft dazu, dass der Graph-Clustering-Algorithmus Simplizes falsch verbindet oder trennt, da er rein topologisch arbeitet.

`SHGO.jl` qualifiziert dies durch einen zweistufigen Prozess:

1. **Topologische Selektion:** Identifikation von Regionen, die *potenziell* ein Minimum enthalten.
2. **Dynamische Verifizierung:** Lokale Optimierung liefert den "echten" Zielpunkt. Nur Punkte, die zum selben Minimum konvergieren, bilden ein Basin.

Dies führt dazu, dass `SHGO.jl` bei Funktionen wie dem **Six-Hump-Camelback** eine weitaus höhere Stabilität bei der Identifikation der exakten Anzahl von Minima (6/6) aufweist.

## 4. Performance & Skalierbarkeit

Durch die Nutzung von Julia ist `SHGO.jl` in der Lage:

* Die Geometrie-Berechnungen (Simplex-Inzidenzen) in Nanosekunden durchzuführen.
* Große Gitter (`n_div > 50`) zu verarbeiten, ohne dass der Overhead der Python-Objektverwaltung den Speicher füllt.
* Gradienten-Informationen via `ForwardDiff` (geplant) direkt in die Topologie-Analyse einzubeziehen.

---

### Fazit

`SHGO.jl` ist nicht nur eine Portierung, sondern eine **Weiterentwicklung**. Es adressiert die konzeptionellen Schwächen der SciPy-Implementierung (Clustering-Logik) und nutzt die modernen Möglichkeiten von Julia, um aus einem Optimierungswerkzeug ein Analysewerkzeug zu machen.