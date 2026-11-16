# N03: separate approximate numerical STEP artifacts

Piece: STEP writer and its numerical export contract. No existing caller migration.
Exact `write_step`, `StepReport`, `StepResult`, `Answer`, `BoundedEvidence`, their
proof requirements and numerical-body storage remain unchanged.

New opt-in Rust entry:

```rust,ignore
catallax_step::numerical::build(&body, &options, &mut decider)
    -> Result<NumericalStep, Refusal>
```

`body` is N01/N02's checked `NumericalBody`. The planar-line family is one retained
closed shell with finite rational line edges, plane frames and explicit outer
and inner loops. Source line/vertex and plane/vertex residuals need not be
exactly zero. The decimal payload must instead pass numerical admission under
the body's unchanged stated tolerance. Trims use inverse coordinates on the
actual written orthonormal frame, rather than equating source UV midpoints with
exact real parameters. No world-space source point or carrier is healed.
Projected finite-line parameters must retain the emitted edge's orientation;
reversed or collapsed walks refuse even when a short edge fits inside tolerance.
Native spline edges and rational patches use the separate
[curved numerical path](2026-11-15-numerical-curved-step.md), including flat
splines retained as native spline surfaces. Both paths retain explicit edge/loop
topology. The original N03 qualification's exact-zero source-gap restriction is
superseded; its raw reference checks remain evidence for their original files.

The private artifact bytes are STEP AP242, explicitly labelled `catallax numerical
approximate` in `FILE_DESCRIPTION` and the solid name. Source name, source local
tolerance, check/repair counts and `complete_body_bound=unestablished` are recorded
in the file. No global uncertainty about an exact source model is emitted.

`NumericalStepReport` retains the full generating definition and stable keys,
N01 source evidence/tolerance/repair history, separate evidence from numerical
admission of the selected decimal payload, explicit scalar-rounding limit and
actual maximum scalar-rounding error. A scalar rounding enclosure is not a
whole-curve, whole-face or whole-body bound. Both source and payload evidence
remain numerical, including at zero residual. Existing generic-constructor trust
limitations from N01 still apply; this producer grants no exact/bounded authority.

On the original planar-line path, oblique `DIRECTION` ratios are scaled by their maximum absolute component before
decimal selection. The new report field `max_direction_ratio_rounding_error`
records their dimensionless rounding error separately from written lengths;
it is neither an angular nor a spatial bound. STEP normalizes direction ratios.
Oblique lines use unit vector magnitude and retain their explicit finite edge
vertices and orientation, with a corresponding change of local parameter scale.
Coordinate lines keep the existing magnitude representation.

`NumericalStepOptions` fixes the product name, 1..17 significant digits and a
positive scalar-rounding limit before export. Defaults are 17 digits and 1e-7 mm
per written length scalar, independently of source tolerance. No retry increases
these values or any global cap. One active request allowance covers serialization,
payload admission and internal byte checks; without an active scope N01's
unchanged 20M default applies. Caller precision limits remain effective through
the payload's N01/N02 admission.

Decimal payload checking retains the same vertex/edge/loop IDs and derives
finite parameters and UVs on the written planar frames. It rejects collapsed
edges and malformed/crossing/touching hole loops after rounding, rather than
using source checks as file checks. This is a checked decimal-value model, not a
new global self-intersection certificate. The byte buffer also meets the existing
planar seam syntax contract. No curved pcurve-planning certificate is reused.
The oblique payload reconstructs STEP's normalized directions and projects its
written reference direction perpendicular to its written plane normal, following
OCCT's `StepToGeom::MakeAxis2Placement` and `gp_Ax2` construction. It independently
rechecks that selected decimal geometry through N01/N02 at the fixed tolerance.

Callers obtain bytes via `artifact.bytes()` and provenance via `artifact.report()`.
After writing, `artifact.check_written(&reread_bytes, &mut decider)` checks complete
byte identity to that checked artifact and rereads its seam syntax. Building an
artifact does not itself claim a disk write or OCCT acceptance. An altered,
truncated or unrelated file refuses. External I/O failures remain the caller's
I/O errors; no output path is mutated by a failed artifact build.

The typed refusals use `step.numerical-options`, `-unsupported`, `-gap`,
`-rounding-limit`, `-rounding-topology`, `-file-number` and `-file-changed`, plus
existing N01/N02 topology, trim, sample, precision and work refusals. No partial
artifact is returned. Independent raw reference import/edit validation is separate
and applies only to the exact file hashes it checked.

`step.numerical-gap` still applies when exact-model embedding cannot recover
finite source parameters. It no longer requires exact-zero source incidence
when serializing an already checked numerical body.

[The N03 report](../phase1/N03/REPORT.md) records the public controls, actual byte
reread, raw OCCT check and subsequent cut, with remaining scope limits.
