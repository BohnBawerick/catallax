## 2026-11-03: exact equal-cylinder crossing bands

Pieces: 2, 3, 4, 5, 6, 7; `catallax-features`, `catallax-field`, `catallax-mesh`

Change: `RollingSaddleDefinition::crossing` also constructs `SaddleRing::Crossing { narrow_positive }` for equal-radius cylinders with intersecting, nonparallel axes and equal signed offsets. A complete four-arc selection builds shared pole fans rather than refusing the collapsed contacts. Each band has an exact rational degree `(6,2)` carrier, two Bezier spans, and zero construction deviation. `RollingSaddle::is_edge` and `centres_on` recognize its source and centre ellipses; other saddle families return false.

The chart remains `[-1,3] × [0,1]`, but its two `u` ends are **different collapsed poles**, not a periodic seam. `saddle_rows` recognizes one two-edge loop with shared pole vertices as well as the existing closed-ring topology. Display uses shared pole fans, not zero-area triangles or a welded periodic seam. Exact crossing faces have complete zero-bound junction evidence without a kernel-approximation exception; no positive minimum face width is asserted at a pole.

Contact rails retain the carrier's shared boundary controls. Native queries use an additional exact quadratic representation in the same parameter chart, cancelling positive common factors rather than fitting another edge. Homogeneous cross-product coefficients are proved zero over every complete span before caching the reduction. The new `BSplineCurve::query_span_polynomials` returns those equations with their own positive homogeneous weights. `span_polynomials` and `weight_bound` retain their original stored-control semantics. Exact rigid motion preserves the checked reduction.

`ParameterTable` adds `height_bound` in millimetres. Existing cubic-height tables set it to zero; crossing tables bound both angle and rational height. Their combined cylindrical position error is at most `radius * bound + height_bound`. STEP error accounting includes both terms.

STEP meridian cuts may now pass through an existing spline endpoint when the two incident one-sided control hulls certify a transverse crossing. The exporter retains that topological vertex instead of making a coincident vertex or zero-length split, and assigns each incident edge's seam side from normalized boundary traversal. Compact crossing carriers subdivide only their display proof nets longitudinally; neither the two-span carrier nor the chord target or resource limits change.

Whole-body volume, area and centroid use the exact crossing law and an interior half-angle integration chart that remains regular at the zero-width poles. Integration tolerances, panel limits and work defaults are unchanged.

Why: the two contacts coincide at a crossing pole, so a finite-area spherical corner and the regular ring's positive-width fitting proof do not apply. Rational circular sections and whole-polynomial identities permit an exact reduction instead of an approximation.

Coverage: complete pole-to-pole bands, with a strict centre-ellipse tube-reach proof and whole-control-hull clearance from every other support, which must be planar. Partial pole fans, clipped arcs, mixed ball sides, mixed strip-end constructions in the same step, unproved tube folds, and tangent-shadow continuation remain typed refusals. Full inertia, per-face area, per-edge length and independent whole-body field compilation retain their existing rolling-saddle gates; this slice does not claim those operations.

Migration: exhaustively match the new `SaddleRing` variant; do not treat every non-overlap band as a through ring or identify its parameter ends. Initialize `ParameterTable::height_bound` and include it in positional error accounting. High-level recipe clients need no new operation or flag: retain `catallax-recipe/3` and select all four wall/wall arcs through the existing `between` selector. SolidDraft must update its Catallax dependency, use the existing fillet operation and display options, and keep reporting any optional mesh controls that the kernel marks unavailable.
