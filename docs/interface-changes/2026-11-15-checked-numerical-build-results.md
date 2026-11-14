## 2026-11-15: checked numerical build results

Pieces: numbers, topology, application and CLI verification.

Change: `catallax_app::build_numerical` returns `Result<NumericalBody, Refusal>`.
A successful result retains an ordinary B-rep, its generating definition, stated
millimetre tolerance, completed checks and ordered repair history.

Why: numerical modelling may report checks that establish less than a complete
geometric bound. `BoundedEvidence`, `OperationResult`, exact construction facts,
exact checkpoints and existing exporter interfaces keep their contracts.

Migration: existing callers need no change. New numerical consumers use this
separate result and the definition/check accessors. They must not interpret a
sampled maximum as a bound over a face, an edge interval or a body. There is no
conversion to an exact or bounded answer. The stored model is private so a
numerical result cannot be passed to an existing exact consumer by extracting it.

The definition retains shared vertices, finite edges with parameter endpoints,
and oriented faces with ordered finite u/v loops. Indices reference the retained
tables; nonempty keys are unique within each table. N01 accepts line edges,
planes and embedded rational patches, with exactly one simple polygon loop per
face. General interior and concave trims are supported. Extra loops, unsupported
carriers, invalid indices, unused entities and malformed topology refuse.
Later trim work can extend admission without replacing the definition tables.

The required checks are:

- `topology.structure`: the existing topology checker, excluding only exact
  Incident fact lookup. Numerical admission adds no such facts.
- `trim.simple-oriented-parameter-loop`: rational checks for closure, orientation,
  intersections, overlaps and patch-domain containment of the linear trim.
- `edge.vertex-distance`: the curve's two finite endpoints against their vertices.
- `trim.edge-surface-distance`: nine equally spaced locations per coedge, including
  both endpoints, against the retained linear u/v segment on the surface.

Each distance record retains its curve parameter, surface parameters where
applicable, and a 128-bit outward distance enclosure at that location. The
positive tolerance is fixed before the build; any sampled upper distance above
it refuses. No retry increases tolerance, work or precision. A caller precision
cap below 128 bits refuses. The caller's active operation budget applies; without
one, admission and continuation share a default limit of 20 million work units.
These checks do not establish global self-intersection freedom, separation or
any complete approximation bound.

`replace_surface` retains the finite trims and stable keys, then rechecks the
whole new result. `repair_vertex` measures the actual point movement against the
same tolerance and records the vertex key, before/after positions and movement
upper bound before rebuilding. Both preserve previous repair history. A failed
continuation returns no new body and leaves the old result intact. Parameter
edits are not repairs; no repair is represented by an empty list. Historical
movement bounds are not added into a bound on the current body.

Acceptance examples are executable in the application and wire tests. An
interior trim with sampled maximum 1/4500 mm passes a 0.01 mm tolerance; a
1/45 mm maximum refuses. Zero residuals stay numerical. A missing coedge,
crossing, reversed or absent trim, extra loop, invalid index or exhausted work
refuses independently of other successful checks. A concave six-edge interior
trim retains all six edge uses. Editing a support updates the actual checks;
repairing a displaced vertex records its measured movement.

`numerical_build_summary` and `numerical_evidence_summary` are public Rust wire
adapters. The new `catallax-numerical-result/1` record uses `status: numerical`,
`proof: numerical`, and `exact: false` on success. `proof: numerical` names the
result category; it grants no proof status. Checks retain `kind`, `name`, `scope`,
actual sample locations, pointwise `distance_upper_mm` and `sampled_max_mm`.
Rationals serialize as exact fraction strings. Repairs retain their measured
point movements. Refusals publish neither partial geometry nor acceptance
checks. No existing CLI input format or exact manifest changes.

Numerical STEP export is not implemented here. `step` is always null, including
for zero residuals and repaired results. N03 must add a numerical writer with
its own output checks before declaring an approximate file exportable. N02 and
later numerical builders may rely on the retained definition, tolerance, actual
check records and repair history; they may not infer any unrecorded check.

New refusal codes are `numerical.invalid-definition`,
`numerical.unsupported-carrier`, `numerical.invalid-trim`,
`numerical.unknown-face`, `numerical.unknown-vertex`,
`numerical.repair-tolerance-exceeded`, `numbers.invalid-numerical-evidence` and
`numbers.numerical-tolerance-exceeded`. Structural, precision, work and deadline
failures preserve their existing codes.
