# STEP export blind tests

The blind tests use `docs/pieces/07-step-export.md`, the judge file contract, and the public interfaces of pieces 1 to 5. The loopless-face import regressions additionally use pinned OpenCASCADE through the configured direct translation, which bypasses caller-side processing, sewing and SameParameter; internal repair inside the reader remains unproven.

## Commands

Run the enabled exporter tests, private tests and implementation-independent checks:

```sh
cargo test -p catallax-step --release
python3 crates/catallax-step/tests/verify_exact.py
```

The required CI job discovers all `tools/step-regressions/test_*.py` with
`cadquery-ocp==8.0.1.1.0`. Complete-face and band checks cover spheres, tori,
cylinders, both frustum directions, rational poses, reversed through-hole
walls, enclosed cylindrical voids and algebraic countersinks.
`test_partial_periodic_trims.py` adds a public cylinder/sphere intersection,
touching branch vertices, a closed cylinder/sphere branch, spherical lenses at
identity and a rational pose, and a public block/cylinder union whose two
cylinder band seams each meet two cut vertices. An L-shaped prism minus a ball
at two poses checks that an interior pole moves onto an export-only cap, a
ball cut by a plane passing 1e-15 from both of its frame's poles checks the
retry on the perpendicular frame, one octant of a ball at two poses checks
that a pole on a trim gets another sphere frame, and a dimple touching the outer
boundary of a reversed planar face checks that the split planar loop starts at
its repeated vertex. Each
of these imports the file through two procedures: raw translation and the
standard reader. Each procedure must give a
valid shape with one solid, the expected face, edge and vertex counts, and a
converged adaptive volume within the diagnostic tube of the expected volume.
The raw procedure is the configured direct `StepToTopoDS_Builder`
translation. It bypasses caller-side processing: it never calls `TransferRoots`
or `ProcessShape`, and it keeps SameParameter and max precision off. Internal
repair inside the reader remains unproven. Two scenes in the same file pin the
legacy fallback and check only the export: a transverse ball frame whose pole has no membership
certificate still answers within the request allowance with no added edges, and
touching spherical trims keep a file with no pcurves.
`test_p3_round_charts.py` exports the public P3 filleted boss at both junction
radii and checks raw and standard validity, finite face parameter domains, one
STEP face per source face, converged volume against an independent construction
and three section probes that move with the edit.

`test_step_import.py` reads files through `check-step.py`, which uses only the
standard reader with its default processing. That check requires a valid
shape and the expected number of distinct faces across the whole shape. It
does not count solids or translate raw, and it reports the kernel volume
comparison as unmeasured. The tests also check the regularization report. The parallel-gap cases compare the
standard reader's volume with 168 within 1e-9. On the full slot-cut capture,
whose rounded planar sliver the raw translator cannot resolve, the test pins
the legacy encoding by requiring no added edges.

The Rust `partial` target preserves the 100,000,000-work request allowance,
the legacy encoding of unresolved subprecision trims, and whole-export legacy
topology for mixed conical scenes in both body orders.

Run **all tests**, including the eight cases retained behind individual ignore reasons:

```sh
cargo test -p catallax-step --release --no-fail-fast -- --include-ignored --test-threads=1
```

Seven intersection bodies remain ignored with historical fixture-refusal
reasons; their annotations
are not proof that the boolean engine is absent. The scale-family body remains
ignored because its greedy corner matcher rejects a valid one-to-one assignment
at generated index 17. This publication does not repair these exclusions.

Ignored tests still compile in normal CI. The six independent checks exercise
the STEP reader, de Boor evaluator, rational distance inequalities, quaternion
formula and topology fixtures. An ignored test is not a passing exporter test.
The report distinguishes the original inclusive failures from supplemental
direct intersection and fitting observations.

Every default-budget `export` call also runs the same model, bodies and options
with `max_precision_bits: 64` and `max_work: 0`. This applies to every export
that its test reaches, including each generated family iteration. The original
inclusive run stops the scale family at its first failing assertion and cannot
construct the boolean-dependent intersection fixtures. The zero-work run may
answer within tolerance, answer with the distance mark, or refuse atomically.
Default-budget calls must still answer.
Zero-work answers use the shared text, rounding, topology, vertex-distance and
surface-distance checks, plus 1,000 exact samples per spline span against both
of its model supports. Added seam geometry and ring origins are checked too.

## Coverage

| Contract cases | Test target | Checks |
| --- | --- | --- |
| 1 | `contract` | Default options, exact header and fixed entities, all fourteen entity counts, 170 total entities, corners and exact report |
| 2 | Every export via `common::reader` | Consecutive IDs, existing references, depth-first first-visit numbering, no unreachable entities; unselected store entries do not change output |
| 3 | `tools/step-regressions/test_periodic_bands.py` | Complete cylinder/frustum bands, ring and seam pcurves, outward and reversed-hole orientations |
| 6 | `contract` and `tools/step-regressions/test_periodic_bands.py` | Cone apex starts the band seam, no vertex loop, rounded angle; raw import of pointed cones and drill tips |
| 4, 5 | `tools/step-regressions/test_loopless_faces.py` | Full sphere and torus raw import, poles, seams, validity and analytic volume at two poses |
| 7 | `contract` | Apostrophe, backslash, BMP and supplementary Unicode, control characters, empty name, literal escape-looking text and syntax-looking strings |
| 8 | `contract` | Empty list, empty bodies, mixed empty and solid bodies |
| 9 to 11 | `contract` | Ordered bodies/lumps/shells/faces/bounds/coedges, void reversal, touching lumps with shared PointIds but separate STEP points, shared sphere surface and circle |
| 12 | `deviation` | Five-digit thirds require a distance mark; seventeen digits meet the default tolerance |
| 13 | `deviation` | Named quaternion pose plus 1,000 random rational poses and sizes; independently computed corner positions and plane normals |
| 14 | `deviation` | Exact `10^30` dimensions, inexact `10^30/3`, and 1,000 random posed sizes across scales `10^-30` through `10^30` |
| 15 | Every export and `deviation` | Rational vertex errors and incident plane/cylinder/sphere distances; cone and torus distances by proven balls; 1,000 posed/scaled primitives including both frustum directions and added vertices |
| 16 | `intersections` | Near-tangent irrational lines and coaxial irrational circle centres remain analytic; Steinmetz ellipses remain shared analytic curves |
| 17 | `intersections` | Crossing cylinders and offset cylinder/sphere branches; clamped spline structure, endpoint decimal equality, closed flags, per-edge spline ownership, 1,000 exact samples per nonzero span, default and `10^-12` tolerances |
| 18 | `intersections` and every default export's zero-work companion | 64-bit runs with 10,000 or zero work meet tolerance, return a distance mark with a valid bound, or return `step.budget-exhausted` or an unchanged registered piece 1/2/4 refusal code without altering output |
| 19 | `contract`, `intersections` | Repeated exports and eight cloned-model threads compare bytes, reports and marks, including fitted branches |
| 20 | `contract` | Expired deadline returns `numbers.deadline` without altering output, or exactly the normal bytes/report/marks |
| Refusals and other rules | `contract`, `refusals` | `step.io-error` with exact error detail after immediate and partial writes; `intersect.unsupported-curve` propagated unchanged and atomically; zero digits clamp to one; uncertainty rounds upward; body marks are not copied into exporter marks |

Branch-node regressions preserve the existing shared physical vertex; periodic chart cuts may add other vertices and seams. Surface and curve definitions are matched to the source model, with ordering and sharing checked by IDs, not just aggregate counts. Loopless torus seams carry two pcurves on their unchanged analytic circles. A loopless sphere uses a sole vertex loop for the complete natural domain.

Complete cylindrical and conical bands with two vertex-free ring loops retain
their existing seam encoding. A cone face bounded by one such ring and a vertex
loop at its own apex is a band whose seam starts at the apex. Regular partial cylinders and pole-free spherical
faces use exact shared-edge cuts and bounded, parameter-coupled pcurves. Source
loop-corner incidence is preserved at touching vertices; disconnected chart
components receive distinct face/surface identities. No kernel topology is
modified. An interior spherical pole moves onto an export-only cap behind a
latitude ring. A sphere whose pole lies on a trim is written with the first of a
fixed list of exact frames whose poles are certified outside each of its faces,
or inside one where a cap applies. Elliptic cylinders, cone faces, arbitrary
torus trims, poles without a membership certificate, spheres that no listed
frame admits, touching spherical trims, and
unresolved rounded edges anywhere in the export retain their previous encoding. Uncharted planar caps must keep native analytic boundaries;
opaque boundaries reject the whole export's partial cutover. Complete bands
remain on their earlier path. Distinct endpoints must remain separated by the larger of the requested
distance tolerance and `1e-7`, plus their certified rounding errors.

Longitudinal seam directions come from the certified carrier frame, not from
normalizing an algebraic endpoint-height difference. Full-carrier incidence
and coaxial facts are reused when present; eligibility otherwise retains its
exact proofs. The CLI standard-hole allowance regression stays unchanged.

`step.not-solid` and `step.unsupported-geometry` are not constructible through today's public enums. Their future variants are deliberately not faked. Invalid handles and failed topology checks are programming errors, not refusal fixtures.

Case 18 accepts `step.budget-exhausted` and
registered `numbers.*`, `shapes.*` and `intersect.*` codes, including
`numbers.budget-exhausted` and `numbers.precision-limit` promoted by `trace`.
It rejects all other `step.*` codes, unknown codes and any output mutation on
refusal. A valid but over-tolerance spline fit must still answer with its mark.
`refusals.rs` separately checks unchanged upstream refusal details. The architect
owns the corresponding piece-description changes.

## Exact oracles

`verify_exact.py` uses Python's standard-library `fractions.Fraction` only. It independently checks the 170-entity count, the squared rounding-distance bound `1/30000000000`, the large-coordinate error `10^13/3`, a rational alternating-series enclosure of `atan(3/4)`, the two adjacent decimals around `sqrt(2)`, the named quaternion matrix and the intersection equations.

The Rust spline evaluator uses de Boor's recurrence in exact rationals. Each nonzero span is sampled at `lo + (hi-lo) j/999`, for `j=0..999`. Cylinder and sphere sample distances are checked by squaring the two inequalities around the radius. This checks the actual distance against `max_deviation`, not a residual or piece 2's potentially two-times-larger conversion bound. No floating-point number decides whether a test passes.

The reader's independent fixture is the documented empty STEP file, not output captured from an implementation. Active mutation checks reject dangling references, numbering gaps, unreachable entities, wrong walk order, outside-string spaces and CRLF.

The unit regressions `stalled_algebraic_zero_decimal_uses_an_exact_proof`
and `stalled_tiny_positive_decimal_is_not_snapped_to_zero` exercise a squared
cubic algebraic identity whose enclosures retain positive midpoints. The first
requires exactly zero decimal value and error after a piece 1 proof; adding
`1e-200` must not be snapped to zero when 128-bit refinement stalls.

### Void decision

The void-face rule reverses the model's face sense. Piece 3 case 17's inner faces have `Sense::Reversed`, so their exported `ADVANCED_FACE.same_sense` is `.T.`, their `FACE_BOUND.orientation` is `.F.`, and the wrapping `ORIENTED_CLOSED_SHELL.orientation` is `.F.`. The original piece 7 case 10 expectation was a typo; the general rule and mathematical assertions stand.
Complete periodic bands instead keep their `FACE_BOUND.orientation` equal to
`ADVANCED_FACE.same_sense`, including on void shells, because their parameter
rectangle always retains the counter-clockwise chart walk.

## Non-gating measurements

Compile the blind benchmark without running its boolean-dependent section:

```sh
cargo bench -p catallax-step --features implementation --bench export --no-run
```

Run the primitive and face-count measurements. It currently stops when it reaches the unimplemented piece 5:

```sh
cargo bench -p catallax-step --features implementation --bench export
```

`cargo bench` uses the optimized bench profile, so it does not accept `--release`. The `implementation` feature gates only this executable benchmark. Enabled tests no longer require `--include-ignored`.

The CSV output records case, curve tolerance, face count, file bytes, export nanoseconds, spline count, total control points, reported deviation, mark count and outcome. Parsing and counting run after the timed export. Fixtures include all primitive forms, increasing block face counts, all three boolean operations on piece 5 cases 1 to 40, and finite solids with the side intersections of piece 4 cases 17 to 26 at tolerances `10^-6` through `10^-12`. Cases without quartic branches report zero splines rather than inventing a branch. Boolean refusals are printed separately. There is no time, file-size or control-point-count threshold.

The implementation's direct piece 4 trace/fitter benchmark does not require
piece 5 and can measure cases 17 to 26 separately:

```sh
cargo bench -p catallax-step --bench speed
```

These direct measurements are supplemental. They do not turn an unconstructible
blind intersection body into a pass. The comparison records all observed
figures and the reference writer's missing certified fitting capability.
