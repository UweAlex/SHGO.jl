
# README: Grid.jl (Aktualisierte Fassung)

## 1. Übersicht

Das `Grid`-Modul ist der „Vermessungstechniker“ des Gesamtsystems. Es wandelt das mathematische Kontinuum in ein diskretes, adressierbares System um und legt die Basis für die Kuhn-Triangulierung.

## 2. Kern-Konzepte

### A. ID-Hoheit (Determinismus)

Das Gitter definiert eine bijektive Abbildung zwischen einem **kartesischen Index** und einer **linearen ID** (Integer).

* Jede Ecke hat eine eindeutige Adresse.
* Die ID ist der Primärschlüssel für Kuhn und die `TopicalStructure`.

### B. Isotropie-Garantie

Um „Nadel-Simplizes“ zu vermeiden, berechnet das Grid die Zerlegungen () pro Dimension so, dass die Zellen möglichst würfelförmig bleiben. Dies geschieht auf Basis des Parameters `n_div`.

## 3. Die Struktur `GridStructure`

* **Attribute:** `lb`/`ub` (Schranken), `k_steps` (Intervalle), `dims` (Punkte pro Achse), `linear_indices` (Mapping).
* **Methoden:** `get_linear_id`, `get_vertex_pos`, `total_vertices`.

## 4. Der Isotropie-Algorithmus

1. **Volumen-Analyse:** Erfassung der Box-Ausdehnung.
2. **Abstands-Normierung:** Bestimmung des idealen Abstands .
3. **Schritt-Verteilung:** Dimensionen mit großer Ausdehnung erhalten proportional mehr Stützstellen.

## 5. Schnittstellen

| Empfänger | Information | Zweck |
| --- | --- | --- |
| **TopicalStructure** | `get_vertex_pos(id)` | Erzeugt `Vertex` und `work_heap`. |
| **Kuhn** | `get_linear_id(idx)` | Verknüpft Punkte zu Simplizes. |

## 6. Axiome

1. **Unveränderlichkeit:** Gitter-IDs sind statisch.
2. **Lückenlosigkeit:** Jeder Gitterpunkt besitzt eine ID.
3. **Minimalität:** Jede Dimension hat mindestens 2 Ecken.

---

## 7. Die Evolution der Indizierung (Zerkleinerung)

Das Gitter ist die **Initial-Struktur**. Für spätere Erweiterungen gilt:

### A. Das „Verlassen“ des Gitters

Bei lokaler Verfeinerung (Refinement) entstehen neue Punkte, die keinem festen Gitterplatz mehr entsprechen.

* **Gitter-IDs:** Folgen der -dimensionalen Logik.
* **Refinement-IDs:** Werden dynamisch vergeben und besitzen keine kartesische Koordinate im Ur-Gitter.

### B. Die TopicalStructure als Brücke

Nach der Initialisierung übernimmt die `TopicalStructure` die Navigations-Hoheit.

1. **Start:** Kuhn nutzt Gitter-Mathematik für den schnellen Aufbau.
2. **Betrieb:** Die Analyse-Logik nutzt das `star_mapping` (Graphen-Adjazenz). Es ist unerheblich, ob eine ID aus einer Matrix-Rechnung oder einer lokalen Teilung stammt.

### C. Strategie für zukünftiges Refinement

* Das Grid liefert die geometrische Basis (`lb`, `ub`, `delta`).
* Neue Punkte berechnen ihre Position durch Interpolation ihrer Eltern-Punkte im Simplex.
* Die `TopicalStructure` verwaltet diese Punkte als topologisch gleichwertige Bürger, auch ohne festen Gitterplatz.

---

### Warum das jetzt wichtig ist:

Durch diesen **hybriden Ansatz** (Gitter für den schnellen Start, Graph für die flexible Analyse) ist SHGO.jl sowohl performant als auch zukunftssicher für adaptive Gitterverfeinerungen.

---

**Fazit:** Das Dokument ist absolut noch aktuell und bildet die "Verfassung" deines Grids. Es rechtfertigt, warum du dich für `Dicts` statt für `Arrays` im Manager entschieden hast – eine Entscheidung, die sich jetzt beim korrekten Basin-Clustering auszahlt.