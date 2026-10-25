## 2026-10-19: piece 7 deepened, with a fixed file layout

Pieces: 7, cli
Change: no signature change. The description fixes the header, the twenty fixed entities, the entity mapping, the depth-first numbering, the real-number text, the entities added for ring edges and loopless spheres and tori, the writing of voids and branches, and the three parts of `max_deviation`. `write_step` writes nothing on a refusal.
Why: byte-identical output and blind tests need one layout; readers need every face bounded.
Migration: none for the front end beyond reporting `added_vertices` and `added_edges` from the report.

