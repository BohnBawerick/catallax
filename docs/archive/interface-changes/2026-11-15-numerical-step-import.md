Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Numerical STEP source admission

`catallax_import::numerical::import_step` checks an actual STEP source and its
expected SHA-256, then returns `StepImport::Certified` or `StepImport::Numerical`.
The certified variant retains the previously supported import behavior and body
marks. The numerical variant owns separate `NumericalBody` values and the actual
source-repair record; it provides no extraction into an exact `Model`.

The existing `import_step_bodies` entry now commits all source solids atomically.
Failure preserves the caller's complete model, including its allocation state.
It retains its original result category and does not call the numerical fallback.

## Numerical family and evidence

The fallback handles planar faces with finite straight boundaries after a
recorded source-incidence repair. It reuses the recorded producer's source cap,
nearby-feature checks and actual accepted definitions. It does not substitute raw
source geometry when an internal builder repair has not returned its definition.
The requested numerical tolerance must be positive and at most 0.01 mm, and each
recorded source-deviation upper bound must fit that requested tolerance.

Source incidence history remains in `NumericalImport::record`; the numerical
builder's evidence records its own actual checks. These are separate histories.
The conversion retains the accepted planes, finite world-space vertices and
walked loops. Equal endpoint coordinates share a vertex within one source solid;
finite straight edges share their unordered endpoint pair. This is a restricted
straight-edge conversion, not proximity matching between different sources.

Coordinate planes use a canonical unit parameter frame without moving geometry.
Oblique planes retain their accepted frame. Their inverse orthonormal coordinates
use 128-bit enclosure midpoints as rational numerical UV coordinates. The N01/N02
builder then checks the actual resulting endpoint and coedge samples at the fixed
tolerance. Its existing result and sampling contract applies; this conversion
does not turn sampled distances into a bound between samples or over the body.
Zero sampled residual still returns the numerical variant.

The caller's active work allowance spans the initial certified attempt, recorded
source repair, conversion and numerical admission. Without an active request,
the complete query uses the unchanged 20M default. Errors distinguish source
refusals from numerical refusals and return no partial geometry.

The standalone `numerical_import` example attempts source admission, a pose
rebuild and approximate STEP writes with disk rechecks. It records completed
stages and preserves admission evidence when a later export refuses. A failed
workflow is not counted as a complete import/edit/export qualification. The
[numerical writer's contract](2026-11-15-numerical-step-export.md) governs STEP support.

## Approach and verification

The source repair follows the existing recorded import producer after inspection
of OCCT's `ShapeFix_Face` implementation. The planar inverse follows
[`ElSLib::PlaneParameters`](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_9_2/src/ElSLib/ElSLib.cxx),
which transforms a point into the plane's orthonormal frame. Catallax retains
numerical UV values and independently rechecks their world-space residuals.
The [P04 SMLib port map](../phase1/P04/PORT-CONSUMERS.md) was inspected; this
analytic plane conversion needs no tracer or general projection port.

The focused integration controls cover a repaired source, unchanged supported
inputs, a rigidly rotated oblique source, the inclusive source cap, a larger gap,
the nearby-feature limit, a stricter requested tolerance, a true crossing,
multi-solid rollback, exhausted work, and a pose rebuild with an actual
approximate disk STEP write and reread check.
