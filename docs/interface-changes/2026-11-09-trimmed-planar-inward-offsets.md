## Inward offsets of nonconvex trimmed planar solids

Ordinary B-rep `offset` now accepts negative distances on bounded planar polygonal solids whose complete finite nearest-feature arrangement is proved. This includes nonconvex extrusions, polygonal holes and planar faces retained from Boolean operations. The recipe format and distance semantics are unchanged.

The convex, radial, cylindrical and conical dispatch paths retain their existing reductions. A negative offset that leaves the convex planar family now uses the finite face, edge and vertex distance envelopes already used by planar shells. It selects the source material outside every boundary-distance tube. This gives cylindrical boundaries at reentrant edges and spherical boundaries where a finite vertex contributes. It does not replace finite trim distances with support-plane intersections.

New planar patches retain their source face origins, as in the convex inward route. New curved patches use the existing round edge and corner origins. Meeting envelopes undergo the ordinary exact arrangement and sewing checks. A consumed cavity returns an empty body. Exhausted geometric proofs retain their named upstream refusals.

The spatial reentrant-corner regression builds and exports a valid body, while native volume evaluation retains `measure.budget-exhausted` for an uncovered finite curved trim integral. OpenCASCADE read-back checks its independent `44 - 5*pi/3` volume oracle.

Outward nonconvex offsets, curved trimmed-face offsets, spline normal envelopes, and sharp multi-face thickening remain separate work. The independent field compiler still returns its existing coverage refusal for nonconvex planar offsets. SolidDraft UI behavior is outside this change.

Verification uses `catallax-features --test trimmed_planar_offsets` and `tools/step-regressions/test_trimmed_planar_offsets.py`. Native regressions check material, exact boundary membership, volume enclosures against independent analytic formulas, holes, Boolean-trimmed sources, consumed cavities and rigid placement. STEP regressions import the raw files without healing, check validity and face counts, compare analytic volumes and material probes, and replay the exported question with identical STEP bytes.
