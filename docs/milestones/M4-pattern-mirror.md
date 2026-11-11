# M4 pattern and mirror qualification

The kernel gate combines complete target recipes with body copies and saved feature expansion. It rebuilds each copy after the target's parameter edit and compares face origins, strict STEP bytes, field query results and logical work at 1, 2, 4 and all available kernel threads.

Each target, placement, copy kind and thread count runs as a separate test. Every case compares both parameter states against fresh single-thread builds. The 100M per-snapshot allowance and the existing test hang detector stay unchanged.

| Target | Parameter edit |
| --- | --- |
| P1 bracket | Boss radius 8 to 9 mm |
| P2 collar | Chamfer distance 1 to 0.75 mm |

Each target runs through a two-body linear pattern, a two-body circular pattern about an offset nonunit axis, and a reflection in an oblique plane. The circular cases double all target dimensions and the edited parameter. The same placements expand the complete feature graph into ordinary saved operations with explicit instance ids. Both original and edited recipes round-trip through recipe serialization before evaluation.

The gate checks the expected face count, copy-specific origins and complete per-face origin sets across the edit. STEP output must change after the edit, retain the existing `1e-7` mm export bound without qualification marks, and match byte for byte across thread counts. Field checks certify a copied interior point and surviving named boundary points away from edges. Whole-part distance comparison has its own gate.

The sequential circular-hole case expands three cuts on a square plate. Every instance reads the preceding cut result, keeps its own creating-step name, and retains a selectable wall in the final field result. Its volume enclosure must contain `3600 - 16*pi` cubic millimetres.

P2's reflected feature graph also exercises exact axial separation of coaxial chamfer regions. When world boxes and radial ranges overlap, the fallback projects each circular edge, contact and trim vertex onto the shared axis. Parallel circle planes have constant projection, and the ruled affected region lies inside their convex hull. Strictly separated projection ranges prove disjoint regions. Touching and overlapping ranges retain the refusal; no normalized or approximate axis is needed. Existing box and radial proofs run first.

Run the gates with:

```sh
cargo test --release -p catallax-app --test target_pattern_edits -- --nocapture
cargo test --release -p catallax-app --test pattern --test pattern_edit_sequence
cargo test --release -p catallax-features -p catallax-field --test patterns
```

The existing exact-placement tests also cover touching, overlapping and coincident body copies, copy-zero names, reversed mirror loops and algebraic 45-degree body rotations. Saved feature parameters that lack an exact rational spelling retain `pattern.expansion-exact-parameter`; expansion leaves the source recipe unchanged. Invalid axes, normals and copy counts retain their typed refusals.

The delivery replay compares the 477 pinned earlier cases and the complete P1/P2 body-copy exports and name edits through fresh CLI workers. Raw evidence retained privately. Offset, shell, thicken and SolidDraft interaction have separate M4 acceptance gates.
