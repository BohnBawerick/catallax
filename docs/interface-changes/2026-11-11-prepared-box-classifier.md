## 2026-11-11: prepared box classifier (F6)

Pieces: field

Change: `catallax_field::printing` adds `BoxClassifier`, `BoxStorage` and the public `bounds(plan, dec)`. `BoxClassifier::new(plan, dec)` prepares every live leaf once, and `classify(lo, hi, storage, dec)` returns the `BoxState` of a closed rational box. `classify_box` now prepares a `BoxClassifier`, so it no longer refuses with `field.print-coverage`. Profile sweeps that are not convex extrusions, revolutions, chamfer and fillet regions, drafts and closed meshes are now classified. Its remaining refusals are `field.print-box` for an inverted box, and the decider's work and deadline outcomes.

Every leaf kind that `classify_box` covered before keeps the same interval rule and the same answer. Other leaves are settled where every carrier polynomial keeps one strict sign on the box, or where a mesh box meets no triangle box. One exact point query then gives the state. [FIELD.md](../FIELD.md#prepared-box-classifier) gives the method. Answers do not depend on the worker count or the order of boxes.

Why: PLAN M13 step F6 asks for a prepared box classifier with no coverage refusal on P1 to P4, for the cut-cell grid of M16 and for M19. Before this change `classify_box` refused P2 and P3 and prepared the plan again for every box.

Migration: none. Callers that classify many boxes prepare one `BoxClassifier` and give each worker its own `BoxStorage`. `printing::mesh_plan`, `mass::properties` and the other print and mass routes keep their coverage, answers, work and `field.print-coverage` refusals.
