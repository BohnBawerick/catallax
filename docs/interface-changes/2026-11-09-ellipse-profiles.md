## 2026-11-09: ellipse profile segments

Pieces: recipe, features, application, field, measurement, SolidDraft

Change: `catallax_recipe::Segment` adds `Ellipse` and `EllipticArc`, using a centre, first semi-axis vector and positive second-to-first ratio. The application converter accepts corresponding `ellipse` and `elliptic_arc` SolidDraft segments and projects arc endpoints within the source join tolerance onto exact rational ellipse points.

Why: sketch ellipse profiles must construct exact solids rather than replace curves with display polygons.

Measurement: complete ellipse boundaries on planar caps and matching elliptic cylinders use exact polynomial integrals for volume, centroid and inertia. Elliptic arc length and curved-face area retain enclosed numerical integration.

Migration: exhaustive segment matches must handle the new variants. Normal extrusion supports lines and co-oriented, equal-aspect ellipses, including nested holes and partial ellipses. Unsupported frame mixtures return `sketch.ellipse-frame`; circular mixtures return `sketch.ellipse-circular-mix`. Nonpositive axes return `sketch.ellipse-axes`; inconsistent exact arc endpoints return `sketch.ellipse-endpoint`. Oblique extrusion returns `extrude.ellipse-oblique`. Revolution, path sweep, first-hit extrusion limits and field profiles return named operation-specific refusals. Existing line, arc and circle recipes retain their format and behavior.
