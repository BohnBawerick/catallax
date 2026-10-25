## 2026-10-19: piece 3 check gains fault codes and entity-sharing rules

Pieces: 3, 5
Change: `Model::check` reports eleven fault codes: the earlier `topology.loop-not-closed`, `topology.edge-use`, `topology.euler`, `topology.disconnected-shell` and `topology.missing-incidence`, plus new `topology.body-kind`, `topology.empty-part`, `topology.shared-entity`, `topology.edge-ends`, `topology.unbounded-face` and `topology.vertex-loop`. Missing incidence now also covers an edge's end points on its curve and a vertex loop's point on its surface. The Euler count is stated exactly. Bodies may share entities with other bodies; within one body no face, edge or vertex is used by two shells. `owner_of_coedge` gives the first face added that holds the coedge.
Why: the earlier list could not report a vertex loop alone in a face, a loopless plane face or an entity shared between shells, and ND4 needs the last one checked.
Migration: piece 5 records `Incident` facts for every vertex it adds on an edge's curve (piece 2 does this for `Meet` points), so its results pass `check`. Piece 5 may reuse operand faces in a result.

