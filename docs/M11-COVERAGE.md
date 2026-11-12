# M11 spline-fillet coverage

The production builder rounds one convex edge between two complete single-span
rational patch faces. Each end meets one rational plane. Extra trims,
multi-chart supports, arbitrary end faces and junctions remain outside this
construction. The blend port is partial, not complete.

## What main already proved

- [Body construction](https://github.com/BohnBawerick/catallax/pull/521) rounds
  the M27 ruled loft at radius 0.1 mm in 10,349,665 of 20,000,000 work units.
  Its 18-span face has a two-sided rolling-ball deviation bound of 59 nm.
- [STEP export](https://github.com/BohnBawerick/catallax/pull/523) keeps all chart
  curves, bounds both cap cuts by 0.0375 nm and opens as a valid solid in
  OpenCASCADE without repair.
- [Imported construction](https://github.com/BohnBawerick/catallax/pull/526)
  rounds the same edge on the OpenCASCADE-written loft in 9,158,476 units.
  Its 18-span face has a 25 nm bound. Both bodies keep every face and edge name
  when the radius changes to 3/32 mm.
- The shape-level `m11_fillet_faces` test also proves a face between two
  biquadratic supports. The term `u(1-u)v(1-v)/5` bends their interiors in
  opposite directions. That face uses 22 spans, reaches 80 nm and costs
  16,433,004 units. This is a face proof, not a completed non-ruled body.

The target remains 100 nm, the hard cap 1 micrometre and the application
allowance 20,000,000 units. This change alters none of them.

## Point queries

`classify_point` and `face_contains` use the represented body, including its
physical trims. A trimmed imported patch uses its proved injective projection
and boundary winding. The loft's one-chart atlas uses that chart's embedding
projection rather than requiring affine source parameters.

The cap boundary trace encloses the exact plane-cut graph `tau(theta)`.
It does not substitute the STEP approximation. The existing monotonicity proof
brackets every root for the requested interval of `theta`; interval signs only
discard regions proved empty. Twenty-four bisections tighten each bracket.
An unresolved sign retains the bracket. A query too close for those boxes to
separate retains a refusal.

A blend face needs a checked band domain. The query verifies two exact
restricted carrier rails, two ordered disjoint cap cuts, their endpoints,
physical vertex cycle and orientation. It then tests the cap planes. A surface
match or an origin name cannot grant membership to a missing or reversed loop.

Cap-cut traces currently supply position enclosures and exact endpoints.
Interior exact-point evaluation refuses with
`intersect.fillet-cut-exact-point`; derivative enclosures refuse with
`intersect.fillet-cut-derivative`. A position box is never used as a derivative
bound. Exact incidence in the interior of the carrier's enclosure width also
remains unproved by the existing surface query.

The feature regression builds both public recipes and exercises the native
queries at 1, 2 and 4 workers. It compares complete answers and deterministic
work. Separate tests check cap-cut enclosures against independent height
bisection and reject incomplete, duplicated and reversed band boundaries.

The sequential query matrix records these work counts, independently of the
body-build allowance. Each query keeps its own 20,000,000-unit allowance.

| Point in mm | Answer | Authored loft work | Imported loft work |
| --- | --- | ---: | ---: |
| `(0,0,1)` | Inside | 237,471 | 260,688 |
| `(4,0,1)` | Outside | 512,444 | 437,648 |
| `(1,-1/2,1)`, the removed edge | Outside | 604,906 | 497,391 |
| `(0,0,0)` | Exact boundary | 214,305 | 251,585 |

All seven faces also answer membership of an exact rail point. Exactly the
trimmed support face and the blend contain it. An unresolved interior support
query reports its retained `tau` and `theta` preimage intervals and precision
in bits with `numbers.geometric-gap-unproven`.

## A non-ruled edge

`smooth_recipe` in the feature regression adds a third profile at height 1 mm
to the M27 loft and uses natural cubic interpolation. Its selected edge bends
through `(1,-1,0)`, `(5/4,-1/2,1)` and `(1,0,2)`. Each support has two charts.

At radius 0.1 mm this recipe refuses in the `round` step after 28,255 of
20,000,000 units with `fillet.spline-multichart`:

```text
support chart counts [2, 2]; cross-chart offset-spine tracing and rail clipping
are not implemented; 0 fitted spans, deviation unmeasured; maximum unknown
```

No offset spine or fitted body is published. There is no deviation number to
claim for this case. The missing proof must follow the spine and rails across
the chart join before the existing fit and end-cut construction can apply.
The test compares the whole refusal and work at 1, 2 and 4 workers.

## Blend-port verdict

The [third-party inventory](THIRD-PARTY.md) records these adapted SMLib units.
They generate candidates; production `SplineFillet` still uses Catallax's own
offset-pair proof and Hermite rational-arc fitter.

| Ported part | Gate | Coverage still missing |
| --- | --- | --- |
| `SmConstantRadiusFS` equations and `SmFilletSolver` evaluation, parameter tangents, plane solve and stepping | `catallax-trace --test rolling_ball` | Caller supplies topology and initial second-surface values. Evaluator failures refuse rather than accepting SMLib's failed first guess. |
| `SmFilletIntersector` values, delta values, law point, span accuracy and next point; `SmCrvOnSurf::Evaluate` to one derivative | `catallax-trace --test fillet_intersector` | Offset-domain extension remains `StepRefusal::SurfaceExtension`. Twelve saved cases pin that refusal. |
| `SmSurfaceIntersector::DoPointIntersection` through the fillet stepper, including clipped steps re-solved by the local offset solvers | `catallax-trace --test fillet_trace` | `FlushCurve`, which constructs the fillet surface and rails, is omitted. Seventy-seven saved traces pin the domain-extension refusal. |

The pinned C++ gates compare binary64 outputs bit for bit, including tangent,
touching and singular cases. They are differential candidate gates, not
geometric certificates. Passing them does not complete the blend port or prove
arbitrary spline-face blends.

Catallax's own fitter covers the spline blends measured so far, including the
biquadratic control. General surface fitting is therefore not justified by the
current failures. The next missing work is cross-chart spine and rail proof,
general end trimming and blend construction from the ported trace. Revisit a
general fitting port when a certified spine and rails produce a blend the
current fitter cannot enclose within the unchanged cap.

## Reproduce

```sh
cargo test --release -p catallax-features --test m11_spline_fillet -- --test-threads=1
cargo test --release -p catallax-trace --test rolling_ball --test fillet_intersector --test fillet_trace
cargo test --release -p catallax-shapes --test spline_fillet m11_fillet_faces -- --ignored --nocapture
```

Raw evidence retained privately. Saved replay inputs and approval lists are
unchanged; release replay remains a separate gate.
