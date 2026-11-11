## 2026-11-12: recorded planar source-incidence repair

Pieces: STEP import, feature evaluation, application import records.

Change: `propose`, `apply` and `rebuild` support bounded repair of inconsistent
straight planar source faces. Exact predicates and the unrepaired
`import_step_bodies` entry retain their existing rules.

Why: the [R8 audit](../milestones/R8-planar-self-intersection.md) found that all
745 saved planar refusal pairs contain vertices outside their declared planes.
Such source data does not define one consistent exact face.

### Construction and proof scope

The first candidate refits the affected planes through unchanged vertices.
If that cannot be certified, the second projects each shared vertex onto the
exact intersection of its incident planes. Dependent planes must be consistent.
Vertices used by curved boundaries cannot move. Distinct source coordinates
cannot become one vertex, edges cannot collapse, and shared edge carriers must
remain identical. Straight carriers are refitted through their endpoints when
needed. Candidate selection follows this fixed order, without floating-point
comparisons.

The source interpretation is explicit: `projected_vertex_boundary/1` is the
finite polygon obtained by orthogonally projecting the original boundary
vertices onto the original declared plane. It is neither a recovered intended
design nor a statement of file-supplied uncertainty. Original bytes and their
hash remain unchanged.

For every changed face, complete strictly convex, single-bound source and
result covers prove the two-sided geometric Hausdorff distance. Corresponding
affine maps also bound every original vertex-fan triangle and finite original
line segment against the result. The line-segment endpoints are the orthogonal
projections of the original endpoint vertices onto their original line carrier.
The local upper bound is the maximum of these directly measured bounds.
Intermediate errors are never added or clamped. Unsupported covers, curved
dependencies and missing neighbouring feature proofs refuse.

Every bound must fit both the inclusive 0.01 mm import cap and a tenth of the
source/result neighbouring feature bound. Explicit seam approval cannot raise
the incidence-repair cap. Final bounds are recomputed against the immutable
source after all repairs. Closure, self-intersection, orientation, positive
volume and topology checks still precede publication.

### Wire changes

`GapKind` gains `planar_incidence`. Incidence proposal ids are `p0`, `p1`, and so
on, with the existing hash-keyed `question` spelling. Such proposals are
automatic and accept a decline only. Declining any linked repair leaves that
whole solid unchanged and unbuilt. Existing seam-only decline behaviour stays
unchanged.

Sources needing incidence repair use `catallax-import-proposals/2` and
`catallax-import-repair/3`. The record adds `source_incidence`, with one entry
per changed face:

- `solid`, `face`, `method` (`refit_planes` or `snap_shared_vertices`) and
  `source_definition`.
- `source_discrepancy_lower_mm` and `source_discrepancy_upper_mm`: the measured
  disagreement of raw vertices and finite line endpoints with the defined
  projected source patch. This is separate from geometric deviation.
- `geometric_deviation_lower_mm` and `geometric_deviation_upper_mm`: the
  two-sided distance between the defined source patch and the result.
- `source_to_result_upper_mm`: the direct bound over all retained source
  constraints, also used as the local `face_tolerances` upper endpoint.
- `vertex_movement_upper_mm`: maximum original-to-result vertex displacement.

Distances are outward decimal enclosures in millimetres. In `/3`, the
`face_tolerances` lower endpoint witnesses geometric deviation and its upper
endpoint covers all retained constraints, even when geometric deviation or
vertex movement is zero. The approved limit remains 0.01 mm. The exact rational
feature-size-squared field retains its existing units and meaning.
Change-record coordinate bounds use exact decimal strings where possible and
exact `p/q` strings otherwise. Projection can produce rational coordinates
without a terminating decimal; these are retained without rounding.

`ImportRefusal` gains `planar_repair_not_verified`, carrying `solid`, `face`,
`stage`, `measurement`, optional measured deviation endpoints, `limit_mm`,
optional `feature_size_squared_mm2`, and `reason`. A point/plane residual is
labelled `vertex_plane_residual`; it never stands in for a whole-face bound.
Resource exhaustion remains a separate refusal.
Recipe evaluation maps this refusal to `import.planar-repair-not-verified`.

### Migration

Update exhaustive `GapKind` and `ImportRefusal` matches. Preserve `/3` source
interpretations and all bound fields when saving records; show the source
discrepancy separately from geometric deviation and movement. Call `rebuild`
on reopen. It recomputes every bound and derivative hash before mutating the
destination model.

Clean and seam-only imports still emit the existing proposal `/1` and repair
`/2` formats. Their JSON omits the empty new field. Existing `/2` records remain
valid and are rechecked under their existing source definition. Do not upgrade
a record by changing its tag or inventing source-bound fields. SolidDraft's
rendered viewport and persistence rollout are not qualified by this slice.
