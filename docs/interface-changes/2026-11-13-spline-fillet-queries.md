# Spline-fillet point queries and coverage refusals

The existing native `classify_point` and `face_contains` queries now handle
the physical trims of the authored and imported M11 loft fillets.
Construction, naming, STEP and recipe formats are unchanged.

`SplineFillet::cut_trace` returns an opaque `FilletCutTrace` only after proving
the plane's unique transverse cut. It exposes the checked cut, exact rail
endpoints and interval position enclosures. `catallax_intersect::trace` uses
those enclosures for a cap cut joining both rail endpoints. Unsupported partial
arcs refuse with `intersect.fillet-cut-arc`.

`TracePiece::evaluate_exact` supports these traces' endpoints. Interior exact
evaluation refuses with `intersect.fillet-cut-exact-point`; derivatives and
Taylor orders above zero refuse with `intersect.fillet-cut-derivative`.
The binary64 Taylor filter returns `None` so callers use the checked path.

`boolean.spline-fillet-domain` identifies missing or inconsistent physical
band boundaries. `fillet.spline-multichart` reports each support's chart count,
the missing cross-chart offset-spine and rail-clipping proof, zero fitted spans
and an unmeasured deviation. It replaces the generic elementary-topology
refusal on that unsupported family.

See [M11 coverage](../M11-COVERAGE.md) for the proof scope, limits and the
partial blend-port verdict.
