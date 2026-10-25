## 2026-10-19: piece 3 fixes its builders' geometry and refusal order

Pieces: 3, 5, 6, 7
Change: the piece 3 description now gives each builder's exact surfaces, curves, points, face order, loop cycles and coedge senses in the local frame, the cone surface's apex, axis sign and slope for both frustum directions, the facts each builder must record (with origin names `topology.block` and so on), and how many entities each adds. Builders check `topology.non-positive-size` before `topology.cone-is-cylinder` and `topology.torus-self-intersecting`, and a refused builder adds nothing to the model, the store or the fact book. No signature changes.
Why: blind tests compare builder output field by field, and later pieces read the builders' seams, senses and facts.
Migration: none for signatures. Piece 5 can rely on the listed facts for every operand edge and vertex.

