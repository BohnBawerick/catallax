## 2026-10-26: valid STEP natural boundaries for loopless faces

Pieces: 7, STEP consumers

Change: A loopless sphere now writes one `VERTEX_LOOP` at its south pole.
The sole bound denotes the complete spherical domain. OpenCASCADE's raw
translation constructs one seam and two degenerate pole edges, without healing.
`StepReport` counts explicit file entities, so this face contributes one face,
zero edges and one vertex, with `added_edges = 0` and `added_vertices = 1`.
A loopless torus retains its two analytic circles and one vertex, but writes
`SEAM_CURVE` geometry with paired pcurves and an outward-oriented closed
parameter rectangle. Its counters remain two added edges and one added vertex.
The certified decimal rounding of the torus period contributes to
`max_deviation` through the surface-map bound `major + minor`.

Why: The old sphere seam omitted the pole boundaries. The old torus boundary
left seam selection ambiguous. Both failed independent raw OpenCASCADE import.
These encodings replace the loopless-face rules in piece 7 and make a narrow
exception to ND5's prohibition of pcurves, only for full torus seams.

Migration: Update expected STEP bytes only for outputs containing loopless
sphere or torus faces. Other exports retain their bytes. Do not interpret
`StepReport.edges` or `vertices` as counts of topology reconstructed by an
importer. Model topology and exact geometry are unchanged. This correction
introduces no displaced geometry or gap regularization, so existing
`StepRegularization` records and sidecars remain unchanged. Added entities use
the existing counters and rounded geometry uses the existing deviation policy.
The pinned import regressions bypass automatic ShapeFix processing and keep
SameParameter and tolerance limiting disabled.
