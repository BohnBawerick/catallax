## 2026-10-22: exact conic dimensions and curve closure

Pieces: 4, 5, 7, cli
Change: `CurveForm::Circle.radius_sq` and `CurveForm::Ellipse.major_sq` and `minor_sq` are `Real` rather than `Rational`. `CurveForm` adds oriented `Parabola` and single-branch `Hyperbola` forms. Piece 4 exports `curve_is_closed(geom, curve, dec) -> Result<bool, Refusal>` and boolean callers use it to distinguish ordered open curves from closed loops.
Why: even a horizontal torus section can have irrational squared circle radii. Cone sections also need open conic forms, so curve kind cannot stand in for closure.
Migration: construct exact dimensions with `Real::from` and compare them through `Decider`, retaining exact rational arithmetic when both expressions expose rational values. Primitive geometry definitions remain rational. Analytic storage requires rational fields; otherwise retain the intersection definition. STEP retains the existing rational square-root path and accepts exact real conic dimensions. Parabola and hyperbola forms specify the next pair implementations; this change adds no supported surface pair. The CLI schema and existing curve-kind names remain unchanged. Its `intersection` edge tag still means no recognised conic form; exposing open conics requires a separate versioned exchange decision.

