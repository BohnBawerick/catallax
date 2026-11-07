## 2026-11-07: Exact STEP spline source data and unproved preview

Pieces: 2, 8, mesh

Change: add exact spline source definitions, finite-domain knot extraction and a source-file `PreviewMesh` route. Add exact rational-box position evaluation to `BSplineSurface`.

Why: imported spline source geometry needs an evaluator and display without pretending it has the proofs required for solid publication.

Migration: exhaustive matches on imported `Curve` and `Surface` must handle `BSpline`. Existing body mesh callers retain the default picking types. Call `preview_step` for unverified imported source data, and do not pass its output to certified export or print paths.

`catallax-import::read_step` and `read_step_with_decider` recognize spline curves and surfaces in `MANIFOLD_SOLID_BREP` source files. This extends source geometry recognition, not certified recipe-solid coverage or curved repair.

### Exact source data

Simple `B_SPLINE_CURVE_WITH_KNOTS` and `B_SPLINE_SURFACE_WITH_KNOTS`, and their complex instances with `RATIONAL_B_SPLINE_CURVE` or `RATIONAL_B_SPLINE_SURFACE`, are recognized. Bezier, uniform and quasi-uniform curve/surface encodings generate their exact knot vectors. Decimal controls and knots remain `Rational`; lengths alone are scaled to millimetres. Weights and parameters are dimensionless.

`Curve::BSpline(StepSplineCurve)` and `Surface::BSpline(StepSplineSurface)` retain separate Euclidean controls and positive weights. Surface storage is u-major, v-minor. Closure and self-intersection logical declarations are retained, including `.U.`, but do not authorize geometric decisions. Expanded knot dimensions, monotonicity, multiplicity, active-domain width, rectangular nets and positive weights are checked. Invalid source data refuse as `Malformed`; exhausted caller work refuses as `ResourceLimit`.

The shared `catallax_shapes::spline_tools` module restricts a nonclamped source to its exact active interval with homogeneous endpoint insertion. Tensor restriction applies the same exact operation in each direction. Already clamped data reuse their source storage. This is a change of basis, not a fitted approximation or binary64 round trip. The original file remains immutable and hash-addressed by the existing import protocol. Canonical derivative hashing includes the spline definitions and logical declarations.

The data can construct Catallax's existing `BSplineCurve` and `BSplineSurface` interval evaluators. Those evaluators do not prove regularity, trimming, shared incidence or solid validity. A new `BSplineSurface::evaluate_box` evaluates position over an exact rational box, including non-dyadic finite endpoints, without calculating unrequested derivatives. Its domain and precision checks are shared with `jet_box`; an externally supplied out-of-domain box still refuses.

`build_file`, `build_solid`, repair and recipe publication do not promote general spline data to checked `SplineSupport`. Spline solid construction refuses as `Unsupported` with the source face or edge entity. General spline intersections, orientation and whole-result proofs remain prerequisites for that separate extension. Reading a source file is not a claim that its declared `CLOSED_SHELL` bounds a valid solid.

### Unproved preview

`catallax_mesh::preview_step(bytes, options, dec)` reads the source and returns `StepPreview` with its SHA-256, exact `StepFile` and `PreviewMesh<String, String>`. Picking keys are the original STEP face and edge entity names. Callers retain unchanged source bytes through the existing import-file mechanism.

`PreviewMesh`, `PreviewFace`, `PreviewEdge` and `Triangle` accept generic picking keys. Their default parameters remain B-rep `FaceId` and `EdgeId`, so existing body-preview and certified-mesh consumers keep their current types. No fictitious B-rep handles are allocated for unverified imported data.

The source preview samples each `EDGE_CURVE` once and reuses its vertex indices on every incident face. It uses the existing analytic charts and exact spline evaluators, numerical spline trim projection, constrained parameter-space triangulation and rectangle/triangle subdivision. Holes, reversed bounds and face sense are carried into the triangles. `PreviewOptions` supplies visual targets and a vertex limit; the caller's decider and active operation budget are never replaced or enlarged. Failure returns no partial preview.

Projection convergence, flatness samples, winding and triangles are unproved. There is no conversion to `Mesh`, `Answer<Mesh>`, a solid certificate, checked STEP or print output. Faces without a finite usable trim, unresolved projections, periodic trims requiring an explicit seam chart, collapsed/crossing constraints, exhausted work and vertex/refinement limits refuse instead of showing an untrimmed carrier or dropping a face. A successful preview makes no global accuracy, manifold, material or export claim.

### Exercised paths

The public `tests/fixtures/spline-roof.step` source has five analytic faces, one weighted spline face and one weighted spline edge, with non-dyadic parameter endpoints. The API smoke reads the file, retains its source hash and previews all six faces. The roof centre evaluates to `[0.5, 0.5, 1.1]` mm, and every triangle edge of that fixture has exactly two uses. These fixture observations are not a general certificate.

Regression coverage includes exact nonclamped weighted extraction, a tiny nonzero knot span, U/V net order, source units, malformed data, resource refusals, shared source picking, deterministic output and a trimmed hole with reversed bound orientation. Existing analytic import, certified meshing and STEP export tests remain in their normal suites.

The attributed donor units and licence boundaries are listed in [THIRD-PARTY.md](../THIRD-PARTY.md). Surface fitting is deferred until a new definition-held family needs it. No existing carrier is refitted and the certified STEP intersection-curve fitter is unchanged.
