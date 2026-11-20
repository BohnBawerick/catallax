# Edge copies in Catallax: what exists, who reads it, what the checker checks

Worktree: Catallax repository, commit `bb2f8375 Checkpoint 2: final Catallax v0.4.0 reconciliation (#699)`.
Read-only study. "Inferred" marks anything not read directly.

## 1. The exact model stores no face-curves at all

- `crates/catallax-topology/src/entity.rs:44-62`: `Edge { curve, sense, start, end }`, `Coedge { edge, sense }`. No pcurve field anywhere on edge or coedge.
- So every exact path (primitives, extrudes, analytic Booleans, exact lofts and sweeps) already has one authoritative copy per edge: the 3D curve. The ACIS/Parasolid "exact edge" model is already what Catallax does for exact bodies.
- Face parameters are derived on demand by consumers:
  - Exact tessellation projects 3D boundary samples into the face chart: `crates/catallax-mesh/src/ruled.rs:231` (`chart.canonical(chart.project(...))`).
  - Exact STEP export derives pcurves at write time as inverses of the rounded 3D geometry: `crates/catallax-step/src/pcurve.rs:1-2` ("Parameter-coupled inverses of rounded STEP geometry, not a second rounding of the exact model"), entities written at `crates/catallax-step/src/certified_spline.rs:480-504`, seam/pcurve presence checked in `crates/catallax-step/src/seams.rs:150-216`.
- The numerical checker (`numerical.rs`) does not run on exact bodies. Exact bodies use `model.check` (`crates/catallax-app/src/lib.rs:482-489`).

Conclusion for exact paths: nothing to skip. The owner's proposal is already the design here.

## 2. Two copies exist only in numerical bodies

`crates/catallax-topology/src/numerical.rs:46-99`:
- `NumericalEdge.curve` + `parameters`: the 3D curve.
- `NumericalTrim.uv: [[Rational;2];2]` per edge use: u/v endpoints. Between them the trim is linear, or projected for plane faces (doc at `numerical.rs:65-68`).
- Optional `NumericalNativeTrim { face, boundary, use_index, edge, curve: BSplineCurve }` per edge use (`numerical/native.rs:13-26`): a full 2D spline in the 3D curve's parameter. The shared parameter is a construction convention; the checker samples both at the same t, it does not check reparametrisation separately (`native.rs:203-235`).

### What the checker checks (numerical.rs)

- `edge.vertex-distance` (numerical.rs:709-711): edge ends vs vertices.
- `trim.edge-surface-distance` for ordinary trims (numerical.rs:855-904): 9 points (`SUBDIVISIONS = 8`, numerical.rs:33). At each, 3D point `p = curve(t)` vs `q = surface(uv)`. `uv` is the linear interpolation of `trim.uv`, except for interior points of a plane face with a circle or spline edge, where `uv = parameters(surface, p)` (numerical.rs:868-870, `analytic.rs:102-105`). For those projected trims the check is in effect "is the 3D curve on the plane", not a copy-vs-copy gap.
- Native trims: endpoints and eighths of each knot span, 3D curve vs surface(native uv curve) at the same t (`native.rs:132-235`, check names `trim.native-endpoint-position`, `trim.edge-surface-distance`, numerical.rs:828-840).
- Loop checks on u/v polygons built from trims (numerical.rs:741-810, `analytic.rs:107`, `native.rs:242-380`).
- Face-pair samples at nine edge parameters (doc numerical.rs:166-168). Not read in detail.

## 3. Per path

### Strict numerical STEP import
- 3D curve: read from file; `recognize.rs:609-610` reads only the 3D curve of SURFACE_CURVE/SEAM_CURVE.
- u/v endpoints: derived by inverting the 3D vertex points onto the surface (`crates/catallax-import/src/numerical.rs:399-460`).
- Interior: linear (patch/analytic) or projected (plane). Derived from the 3D curve plus an assumption that the relation is linear.
- Error: zero for plane projection by definition; zero for analytic carriers where the edge really is linear in u/v (line on plane, generator or circle on cylinder, inferred); otherwise unknown, which is what the 9-point check measures.

### Numerical embedding (exact -> numerical; also used for sweeps, lofts, holes, numerical STEP)
- `crates/catallax-features/src/numerical_embedding.rs:1-13`: "Each traced patch-section or patch-pair edge becomes a fitted cubic, and every spline-face edge use gets a fitted native UV curve". Follows OCCT IntTools_FaceFace (one 3D B-spline plus one 2D B-spline per face).
- `fitted_edge` returns `(BSplineCurve, Vec<BSplineCurve>)` (numerical_embedding.rs:688-696): 3D fit and per-face u/v fits from the same trace samples. Each is its own approximation: independent.
- `native_uv` (numerical_embedding.rs:509-530) fits u/v curves for native-carrier edges.
- `numerical_embedding.rs:30-35`: fit densities double "until the binary64 fit meets a quarter of the tolerance. Admission samples every span at its eighths, so the fit is checked at the same places." The fitter refines at the checker's sample points, so a pass says little about the gap between samples.
- Sweeps: `crates/catallax-app/src/numerical_sweeps.rs:120-130` calls `definition_from_exact_model`. Lofts for STEP: `crates/catallax-cli/src/numerical_loft_step.rs:105-120`, same.

### Numerical fillets
- `crates/catallax-shapes/src/spline_fillet/numerical_carrier.rs:200-243`: 3D rails are rows of the blend surface's control net (same knots, weights): on the blend face the rail is an iso-line, derived with zero error (inferred: I did not read the blend-face trim construction, `blend_native`).
- Support-face pcurves (`parameter_rails`) are separate non-rational cubics from packed fit entries (numerical_carrier.rs:207-237), then endpoint-edited by `constrained_rail` (`fillet_numerical.rs:423`, `:836-848`): independent.
- Edited border edges get new u/v endpoints from source native curves (`uv_at`) or plane inversion (`fillet_numerical.rs:806-808`).
- Body rebuilt with `rebuild_with_sampled_carriers` (fillet_numerical.rs:995-1014), so the checker runs again.

### Numerical Booleans
- `crates/catallax-boolean/src/numerical/partition.rs:118-163, 252-294`: section pieces carry their own u/v curves (`piece.uv()`), stored as native trims; source from traced numerical intersection (`crates/catallax-intersect/src/pair/traced/numerical/`). Independent fits (OCCT style). Rechecked via `build_with_native_trims` (partition.rs:277).

### Tolerant STEP import
- 3D curve from file. File PCURVE kept as `source_pcurve` if it matches the sampled edge at 16 stations (`crates/catallax-import/src/tolerant/admission/fit.rs:21, 928-1145`); a mismatch refuses unless boundary repair is on (fit.rs:1139-1145). Otherwise rebuilt: `projected` (plane, keeps source knots), `reconstructed` (poles.rs:181-182), or affine (`step/uv.rs:1-2, 525-545`).
- Both copies retained; `UvSample` carries per-station residuals (`admission.rs:54-68`). Tolerances per edge/face recorded.
- Checks live in tolerant admission, not numerical.rs.

## 4. Consumers

| Consumer | Reads | Evidence |
|---|---|---|
| Exact mesh | 3D curve, projects to u/v | mesh/src/ruled.rs:231 |
| Exact STEP | 3D curve; pcurves derived at write | step/src/pcurve.rs:1-2, certified_spline.rs:480-504 |
| Exact Booleans, measure, SDF field | exact Model (3D only) | entity.rs:44-62 (no pcurve to read) |
| Numerical STEP writer | 3D curve AND trim.uv (as linear PCURVE) AND native trims (as PCURVE); skips planes without native trim | step/src/numerical/curved.rs:124-188, curved.rs:13-40, analytic.rs:815-867 |
| Numerical STEP writer self-check | rechecks rounded pcurves | analytic.rs:620-656, 898 |
| Numerical fillet | source native trims (u/v at border params) | fillet_numerical.rs:408-420, 806-808 |
| Numerical Boolean partition | trim.uv and native piece u/v | boolean/src/numerical/partition.rs:38, 118-163, 488, 543 |
| Tolerant preview mesh | 3D edge samples for edges; trim u/v samples only as inversion seeds | mesh/src/preview/tolerant.rs:130-140, 330-366, 462 |
| Tolerant STEP re-export | retained source PCURVE or rebuilt one | import/src/tolerant/admission/step/uv.rs:525-545 |
| App display of numerical outcomes | none found | app/src/lib.rs:473-477 (inferred: no mesher for NumericalBody; grep found no catallax-mesh use) |
| Drawings/sections | not found in this study | inferred: exact-only |

## 5. Which checks could go

Could skip (derived, zero error):
- Ordinary linear trims on an analytic face whose parameters come from the same construction (line on plane, axis line or circle on cylinder): u/v relation exact. Inferred; no code marks these as derived, so a skip needs a producer-side "derived" tag.
- Blend-face rail trims if they are iso-lines of the blend surface (inferred).

Must stay:
- Plane projected trims: the u/v copy is derived, but this check is the only one that measures whether the 3D curve lies on the plane. It is the gap that matters.
- Embedding 3D fit and u/v fits: both independent; STEP writer emits both.
- Fillet support pcurves: independent, endpoint-edited; read by later fillets and STEP.
- Numerical Boolean section u/v: independent; read by partition and STEP.
- Tolerant retained file PCURVE: independent (authored elsewhere); re-exported.
- Vertex distance and loop checks.

Possibly redundant (inferred): the face-pair distance at the same nine parameters is bounded by the sum of the two edge-surface distances at those parameters (triangle inequality), if both are evaluated at matching t. Not verified.

## 6. Unsafe to skip without changing storage

- Any native trim (embedding, fillet support, Boolean section): the numerical STEP writer emits it verbatim as a PCURVE (curved.rs:146-155), fillets evaluate it (fillet_numerical.rs:806), Booleans partition on it. Skipping its check ships an unchecked independent copy.
- Non-planar ordinary trims: emitted as linear PCURVEs (curved.rs:157-188). The linear assumption is only established by the 9-point check.
- Tolerant retained source PCURVE: re-exported as is.
- Existing weakness regardless of the change: checker and fitter use the same sample points (numerical_embedding.rs:30-35), so no copy is checked between samples (doc numerical.rs:153-155: "pointwise distance enclosures, not bounds between samples").
