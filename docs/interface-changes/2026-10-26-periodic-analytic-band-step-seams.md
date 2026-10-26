## 2026-10-26: explicit seams for complete cylindrical and conical bands

Pieces: 7, STEP consumers, CLI verification

Change: Complete circular-cylinder and conical-frustum faces bounded by two
vertex-free circular ring loops now export one closed parameter rectangle.
The rectangle uses both rings and one longitudinal seam twice. Ring geometry
is a `SURFACE_CURVE` with an exact iso-height pcurve on each incident supported
band. The longitudinal `SEAM_CURVE` carries pcurves at zero and the rounded full
period. Connected bands choose a common physical meridian, reparametrising their
unchanged analytic carriers rather than changing the model.

For a reversed band, `FACE_BOUND.orientation` follows `ADVANCED_FACE.same_sense`
while the chart walk stays counter-clockwise. Both flags reverse together for
void-shell export. This matches the native OpenCASCADE bored-plate convention;
reversing the chart walk instead gives the wrong seam pcurve assignment.

STEP cone pcurves use axial height, not OpenCASCADE's slant-distance parameter.
Decimal axial-coordinate errors are bounded in 3D using `sqrt(1 + slope^2)`.
Angular rounding uses a radius bound over the complete band. These contributions
compose with the existing rounded curve, surface-frame and cone-angle bounds
in `StepReport.max_deviation`. Each repaired band adds one explicit edge and
retains the two added ring vertices. The source topology, exact measurements,
manifest formats and report fields are unchanged. No geometry is displaced and
no gap regularization is introduced; existing regularization sidecars retain
their contract.

Why: Separate noncontractible ring wires without seam and pcurve associations
fail raw OpenCASCADE import, including ordinary cylinders and cylindrical hole
walls. A geometric seam alone is insufficient.

Migration: Accept changed STEP bytes only for supported complete bands and
update added-edge counters. Exact rational cylinders can now report a nonzero
certified deviation because the written period is rounded. This extends the
pcurve exception in [the loopless-face correction](2026-10-26-loopless-step-natural-boundaries.md).
Unrelated files retain their bytes. The pinned regressions cover outer cylinders,
both frustum directions, rational poses and a reversed through-hole wall without
healing, sewing, SameParameter or tolerance limiting.

Known gap: Arbitrarily trimmed periodic faces, elliptic cylinders and singular
apex cones retain their previous encoding; this change does not make them valid
under raw import. A selected file containing such unsupported faces does not
receive a partial band cutover. General certified UV fitting and chart-crossing
cuts remain separate work. Native apex-cone output is not used as a fallback,
because the pinned native writer's own file fails the raw validity check.
