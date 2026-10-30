## 2026-10-30: certified joint pcurves at touching spherical trims

Pieces: 7, STEP consumers, CLI verification

Change: A pinched spherical chart can split into simple cycles and jointly fit its touching pcurves when exact graph-separation and monotonicity certificates succeed. This replaces the touching-sphere fallback described in [regular partial periodic charts](2026-10-26-partial-periodic-step-trims.md) and completes the 262 follow-up named by [the interior-pole correction](2026-10-30-interior-pole-step-retry.md).

Why: Two boundaries can have fourth-order contact. Independently bounded positional fits, even with matching lower endpoint derivatives, do not preserve their parameter-space separation. The joint fit preserves the source loop incidence, keeps physical vertices and space curves fixed, and constrains the two pcurves to opposite sides of a rational polynomial graph.

The writer uses an interior meridian for this class and separates repeated-vertex walks before assigning chart components. At a touching corner it subdivides the first pcurve span and applies bounded opposite control shifts. The common parameter vertex and the original span endpoint do not move. Exact Bernstein signs prove separation, and monotonic graph coordinates exclude new self-intersections with the remaining spans. Optional proof refinement is finite and observes cancellation. An unresolved certificate keeps the existing legacy fallback; arithmetic-budget and deadline failures still propagate.

The UV error includes the maximum actual control displacement, including decimal rounding. The surface projection metric carries that addition into the STEP deviation bound. Neither the requested tolerance nor the raw importer's settings change. Pending pcurve representations are completed only after the joint fit, before STEP graph serialization.

Migration: No recipe, question, native topology or property interface changes. Native answers retain exactness, marks, properties and topology. The approved `invalid-to-valid-step` correction for 262 uses the existing authenticated base-manifest path; 431's accepted bytes remain unchanged. Only the four already-permitted STEP-report fields may differ. The active known-raw-invalid list becomes empty without changing the historical ruling or its fingerprint. Release-correction files remain release work.

The private property audit first tries its existing cheap Gauss-Kronrod sequence. If that sequence does not converge under the unchanged criterion, it uses the prescribed `1e-10`, `1e-12`, `1e-14` sequence. The converged sequence supplies its own spread; an unconverged coarse spread cannot widen the volume comparison. Raw single-solid validity, deviation and area-times-tolerance checks remain hard requirements.

The corrected public touching-trim regression checks both disjoint source components. Each must be raw-valid and standard-reader-valid, with its own independently derived volume comparison. The single-solid test helper retains its one-solid assertion for other fixtures. No application companion or timing-policy change is included.
