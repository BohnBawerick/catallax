## 2026-11-22: STEP assembly exchange

Pieces: 7, 8

Change: `catallax_step::assembly::write_step_assembly` writes a rooted AP242 product graph with rigid parent-relative placements, product and occurrence names, part colours and the existing semantic PMI. `catallax_import::assembly::read_step_assembly` returns that graph with recognized source solids and part-scoped PMI attached to the original face keys.

Why: item 3.2 needs whole products in one file, with shared part definitions and the placement of each occurrence.

Migration: callers use `Assembly<Part>` for a product graph and keep using `write_step` or `write_step_with_pmi` for individual parts. Build imported solids in definition coordinates and resolve PMI through each part's source face keys. Apply the occurrence poses when using those solids in an assembly. The part APIs keep their current contracts. See [assembly and exchange](../milestones/3.2-assemblies-and-exchange.md).
