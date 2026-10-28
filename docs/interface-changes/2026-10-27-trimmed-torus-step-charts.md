## 2026-10-27: bounded STEP charts for circle-trimmed tori

Pieces: 7, STEP consumers, CLI verification

Change: The regular partial-chart export of
[2026-10-26-partial-periodic-step-trims.md](2026-10-26-partial-periodic-step-trims.md)
now also admits a trimmed torus face whose every trim is an explicit circle. Such
a face is cut along an exact tube-circle meridian and written with bounded
degree-seven pcurves on its `TOROIDAL_SURFACE`, like the regular cylinders and
spheres. Scenes that previously kept the whole export on the legacy encoding
because of such a torus, for example the public P3 filleted boss at both
junction radii, now write raw-valid finite charts for the torus, the boss band
and the four partial block-round cylinders.

The written toroidal chart puts its v seam on the outer equator. A torus face is
accepted only when that whole circle is proven outside the face: no trim arc
meets it, and the meridian's crossings alternate from entering to exiting in
tube-angle order, so the meridian arc through the equator is exterior. A trim
on the outer equator, as on a convex rounded edge, or an opaque torus
intersection trim keeps the whole scene on its legacy encoding and its bytes.
Meridian crossings on a torus are ordered by exact tube-angle comparisons
instead of height. No source topology, carrier, tolerance or request budget changes.

Why: The legacy encoding of P3 had six curved faces with unbounded cylindrical V
domains and a full periodic torus chart. The non-repairing raw translator
rejected them and imported 17578.905 mm^3 instead of 14428.640 mm^3.

Migration: None for callers. `StepReport` fields and manifest formats are
unchanged. On P3 the STEP file has one face per source face; the boss band and
the torus band gain one seam edge each, and the three rings gain one origin
vertex each (`added_edges` 2, `added_vertices` 3). Among public recipe fixtures,
only P3 and the fine-mesh filleted boss change bytes, both from raw-invalid to
raw-valid. Known follow-up: these circle pcurves need many fitted spans, and the
P3 export costs several seconds of fitting in the shared pcurve fitter.
