# Piece 8: STEP import and bounded gap repair (contract)

Status: `catallax-import` implements the reader, body construction, repair policy and change records described below. The file-carrying feature and application build entries evaluate `Operation::ImportStep` for one solid, with saved approvals and declines, retained source-relative face tolerances and immutable source verification. SolidDraft preview, persistence and checked-export integration remain separate companion work. See the [local-tolerance contract](../archive/interface-changes/2026-11-02-planar-import-local-tolerances.md) and its [source-incidence extension](../archive/interface-changes/2026-11-12-planar-source-incidence-repair.md).

## Application roundtrip slice

`crates/catallax-app/tests/import_roundtrip.rs` exports a 2 × 3 × 5 mm block and a radius-3, height-7 mm cylinder with `catallax-step::write_step`, then calls `catallax-import::propose`, `apply` and `rebuild` and `catallax_app::build_incremental_with_files`. These entries already existed; this slice adds application-level acceptance evidence, not another importer or a new wire format.

Both sources need no repairs and produce no derivative. Exact carrier, edge-curve and vertex coincidence queries compare the application body with the source primitives. The checks also require sound topology, the file's face origins, matching exact bounding boxes, overlapping volume enclosures narrower than 0.000001 mm³, and a valid checked STEP re-export. Import movement is zero for these fixtures. The source export's reported bound remains a separate quantity and is not inferred from volume agreement.

Display is checked through the application entry at kernel defaults and at SolidDraft's automatic coarse, standard and fine controls, with chord height `0.002 × body diagonal × {4, 1, 0.25}` mm and a 15-degree surface-normal bound. Every face must meet both requested bounds within the existing application work guard. This tests Catallax's controls, not SolidDraft's f32 conversion, rendered viewport or complete application workflow. No outside STEP fixtures are used.

Still refused or unproven by the clean roundtrip slice: splines, curved gap repair and multiple solids in one recipe import step. Spheres, tori and cones of any semi-angle now build (see the surface list below), but this slice does not exercise them. The lower-level importer can return separate bodies; the recipe application does not merge them. Seam exchange contracts and field construction from imported faces are unchanged.

SolidDraft uses `propose` for the preview and `apply` for the repair record, saves unchanged source bytes and their SHA-256 plus the recipe decisions, and supplies those bytes through `ImportFiles` to `build_incremental_with_files`. It must use each proposal's exact `question` string for saved decisions, not construct its own gap-index spelling. `rebuild` checks the persisted repair record on reopen. Recipe values are `decline` or `accept:<tolerance-mm>:<approver>`. Automatic repairs save no acceptance; user acceptance applies only to asked seam proposals.

## Nonzero planar repair qualification

`crates/catallax-app/tests/import_step.rs` exercises nonzero repairs through proposal, application, serialized recipe and repair record, rebuild, cold application rebuild, retained-checkpoint display and checked STEP output. Cases include a 0.004 mm automatic repair, the inclusive 0.01 mm cap, a just-over-cap saved approval, a 0.02 mm saved approval, saved declines, and the inclusive tenth-of-feature limit with a just-over-limit approval refused.

Seam-only repair records use `catallax-import-repair/2`. Each changed convex planar face has an enclosure of its two-sided Hausdorff distance to the complete source trimmed face. The proof checks straight-edge incidence, planarity, ordered strict convexity and complete polygon coverage. Exact vertex-to-trim distances provide lower witnesses. Distance to a convex polygon is convex, so the same witnesses bound every triangle in a complete fan in both directions. This is not a relabelled vertex displacement or carrier-plane residual. Faces with holes, nonconvex covers or unproved adjacent feature bounds refuse rather than acquire a certificate.

The attained upper bound is the local tolerance. It must fit the approved file limit and a tenth of a conservative nearest-feature bound, including narrow skew faces and source/repaired neighbours. A user can approve a larger file limit, not override the feature limit. Bounds are remeasured against the original source after all repairs, never summed or clamped.

`Built::import_record(step)` retains the source-bound record with its immutable import body and recipe, including checkpoint reuse. Rigid placement preserves distances. These records do not certify later modified faces, a whole solid against an intended design, display rounding or exported output. STEP's existing error report remains relative to the repaired model. To bound output against a source face, a consumer needs the applicable source-face record, output correspondence and the separate output bound.

The standalone persisted-artifact smoke approved a 0.02 mm seam and reread the saved recipe, unchanged source, `/2` repair record and derivative. It produced four exact `[0.02, 0.02]` mm face enclosures, a six-face, twelve-triangle OBJ and a checked STEP re-export with zero output error against the repaired model. The separate approval regression grants 0.03 mm while requiring the attained face bounds to remain 0.02 mm. These are API and artifact checks, not SolidDraft viewport or independent OpenCASCADE validation.

## Inconsistent planar source data

The recorded repair entry also repairs straight planar faces whose boundary vertices do not lie exactly on their declared planes. It first tries refitting planes without moving vertices, then exact shared-vertex projection. Changed faces need complete convex covers and source/result neighbouring feature proofs. Every final bound stays within 0.01 mm and a tenth of the feature size. Curved dependencies, nonconvex or holed covers, inconsistent vertex constraints and over-cap bounds refuse with named stages and measured bounds where available.

These repairs use `catallax-import-repair/3` and an explicit `projected_vertex_boundary/1` source definition. The immutable original vertices and finite line pieces are bounded separately against both that projected patch and the result. Records distinguish source-data discrepancy, two-sided geometric deviation, vertex movement and the direct bound over all source constraints. The [interface extension](../archive/interface-changes/2026-11-12-planar-source-incidence-repair.md) defines those fields and migration. Exact predicates stay unchanged. The shared builder retains the separate bounded spline-edge repair and body marks described below. Both paths reserve limb-based exact-arithmetic work; no allowance is raised.

## Supported first path

Files: AP203/AP214/AP242 `MANIFOLD_SOLID_BREP` with one `CLOSED_SHELL`, and `SHELL_BASED_SURFACE_MODEL` containing open or closed shells, in millimetres, inches or metres (exact conversion).

Surfaces: `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`. Curves: `LINE`, `CIRCLE`, `ELLIPSE`, and `SURFACE_CURVE` or `SEAM_CURVE` read through their 3D curve (parameter-space curves are ignored). A `CONICAL_SURFACE` at exactly 45 degrees in a degree context is built exactly as declared. At any other semi-angle, in degrees or radians, tan(semi-angle) is irrational, so no decimal rim away from the placement plane lies exactly on the declared cone. Such a cone is built through the face's axial rim sections, including partial circles, or through its rim and its apex. The apex is a `VERTEX_LOOP`, or the end of a seam on the axis. The slope per retained axis multiple is rational. A raw axis with an irrational length remains exact; the native physical slope divides by the exact square root of its squared length. The whole-face deviation from the declared cone is certified. Endpoint heights, circle projection bounds and positive spline control hulls bound the finite trim domain. The largest radial miss at either domain end gives the upper bound. For exact coaxial rim witnesses, the lower bound is an end's distance to the declared double cone: the smaller of its misses from the two nappes' lines in the meridian plane, times cos(semi-angle). Rounded noncoaxial rims use a zero lower bound and must pass numerical boundary admission if exact incidence fails. The cone is accepted automatically only when the upper bound fits the 0.01 mm automatic cap and a tenth of the face's nearest feature. The feature bound includes the physical axial span, nonzero rim radii and nonzero trim endpoint chords. Numerical cone continuation validates the native carrier and a finite periodic chart strictly above the apex. Its residual allowance is the cap or tenth-of-feature limit minus the fitted-carrier deviation; both allowances are recorded. The deviation is recorded in `face_tolerances`. Beyond either limit the file is refused as `not_solid`, with both bounds in the detail. A sphere or torus face whose bounds reduce to seams, or to one `VERTEX_LOOP` on the surface, is the complete surface. Circles on spheres and torus parallels and meridians are checked exactly. A `VERTEX_LOOP` elsewhere is refused. Several `MANIFOLD_SOLID_BREP` entities import as separate bodies, which acceptance test 2 needs. Bounded rational spline extrusions and B-spline revolutions read as exact tensor surfaces. Other swept profiles, offsets, voids and assemblies still reach typed refusals where their reader or admission path is missing. A face with a single `FACE_BOUND` and no `FACE_OUTER_BOUND` is read with that bound as its outer boundary, as Catallax's own exporter writes it. Decimals read as exact rationals; no float enters a decision.

Face keys are the `ADVANCED_FACE` entity names, as in `docs/RECIPE.md` section 4.18.

Numerical bodies with plane, cylinder, cone and retained ring-torus faces and line or circle edges export through the analytic STEP writer.
Cone pcurves retain axial height as their STEP v coordinate, following OCCT's `GeomConvert_Units::DegreeToRadian` conversion to generator distance on import.
The writer checks the actual decimal placement and semi-angle against each emitted pcurve and edge.
Semi-angle rounding is recorded in radians separately from length errors, with the maximum distance from both endpoints of the enclosed source angle to its written decimal.
Cyclic ordering of a periodic loop preserves the two seam uses on reversed cone faces.
The output remains marked approximate, with its numerical checks and `complete_body_bound=unestablished` retained.
Retained native B-spline trims on cones refuse with `numerical.native-cone-domain-unverified` until complete positive-height and periodic-domain checks cover those finite curves.

Retained ring tori use their original angular surface definition. Finite charts must stay below one turn in each coordinate, or provide a rectangular single-period cell with explicitly paired opposite seam uses. Holes must lie strictly inside each identified period. Every repeated edge needs opposite uses translated by a complete parameter period; a valid first seam does not authorize other wrapped contacts. More general periodic trim layouts remain unverified. Circular rim pairs with matching retained phase anchors can acquire a recorded `paired_torus_iso_circle/1` seam, without moving the source vertices or curves. The writer rotates only the closed walk to preserve OCCT's two seam-chart selections.
Horn and spindle sources retain `import.toroidal-sheet-unverified` on this higher-authority path; selecting only an outer sheet would discard original geometry. Arbitrary native trims retain `numerical.native-torus-domain-unverified`. These are certificate limits, not diagnoses that the source surfaces are invalid.

The separate tolerant native projector retains the source surface and original parameter units. Authored PCURVE hints select the containing source span at each sampled station. An unusable optional hint, including one that cannot compile, preserves the previous fit; work and deadline refusals remain hard failures. Required 3D curve and surface compilation still propagates its own refusals. Cached torus and conic-extrusion guesses initialize projection when needed, including spindle inner-sheet and signed-offset alternatives. Admission still uses the actual evaluated sample, uncertainty, residual and convergence record, not the guess or a closest-distance claim.

Native projection tries the closed-form inverse on unwrapped planes, cylinders, cones, spheres and tori after checking the source seed. Direct samples share the angular enclosures across coordinates. A candidate must still meet the original enclosing residual test. Closed trims on periodic surfaces keep both iterative seam anchors: snapping a reconstructed full-turn cone rim to exactly one period can make the default STEP transfer reject its shared edge. Ambiguous, singular or off-surface projections also retain the parameter grid, Newton updates and line search. Finite chart limits do not need an origin sample or Jacobian to define the search envelope. A bit-identical line-search location reuses the best sample while still checking the deadline.

Sampled projection, fit, sewing and STEP trim reconstruction use the same Rust math routines for norms and planar trace trigonometry on Linux and macOS. A one-bit difference in a residual can change which Newton step is retained, affecting work and exported parameters. The platform-independent residual regression pins these bits without changing the exact enclosures or promoting sampled results to certificates.

For an offset of a regular ring-torus basis, the native source equation keeps the signed minor coefficient `r + d` in the original UV. It does not construct an equivalent checked torus, fit a replacement surface, reverse the chart or acquire an embedding certificate. Horn/spindle bases retain their oriented derivative-normal expression. These equations do not change tolerance limits or the existing RationalPatch proof.

Tolerant STEP export keeps `FILE_DESCRIPTION` short and stores the complete evidence in ordered `CATALLAX_REPORT_JSON` comments. Concatenating their payloads reconstructs the JSON report, including per-entity tolerance, repair and check records; JSON slash escaping prevents embedded names from opening or closing STEP comments.

### Multiple numerical face bounds

The retained analytic admission path accepts multiple interior face boundaries
through the shared numerical builder. Holes must be disjoint, strictly interior
and not nested. Complete planar circular outer rims and holes receive continuous
projected-disk clearance checks under the stated tolerance; general curved
parameter walks retain their sampled checks.
The [numerical trim contract](../archive/interface-changes/2026-11-14-numerical-trim-sewing.md)
defines the scope, padding and evidence names.

### Surface bodies

`StepFile.surface_models` retains surface models separately from `StepFile.solids`.
Each surface model builds as `BodyKind::Surface`, with each shell in its own lump and no assertion of enclosed material.
A closed shell declared inside a surface model stays a surface body.
Boundary edges have one use; interior edges have two opposite uses.
Surface validation checks boundary cycles, edge connectivity, vertex fans and the Euler characteristic with boundary components included.
The display mesh must match the B-rep boundary segments, so permitting an open sheet does not permit missing display facets.

Surface bodies support area, bounds, edge queries, display, picking and STEP export through the ordinary entries.
Volume, solid centroid, inertia and mass properties refuse with `measure.not-solid`.
Boolean operations refuse with `boolean.not-solid` before identity, empty or separated-operand shortcuts.
Solid-field compilation refuses with `field.import-domain`.
STEP output uses `SHELL_BASED_SURFACE_MODEL` with open or closed shell entities determined from actual edge use, and reports zero solids for a surface-only body.
The command manifest reports `result.kind: "surface"`.

The lower-level file entry keeps mixed solid and surface inputs as separate typed bodies.
A recipe that collects both kinds returns a surface boundary collection without asserting combined material; selecting a source solid still returns that solid.
Surface models do not enter solid seam-closing or binary-noise repair.
The numerical fallback refuses a continuation that would omit a surface model rather than return only its sibling solids.
The geometry admission and tolerance requirements for each face remain in force.
Acceptance entries are `catallax-import --test surfaces`, `catallax-app --test surface_import`, and the surface command case in `catallax-cli --test verification`.

### Spline source recognition and preview

`read_step` and `read_step_with_decider` additionally read simple and complex knot-bearing B-spline curves and surfaces, positive rational weights, and Bezier, uniform and quasi-uniform encodings.
Exact rational knot extraction clamps the active domain without fitting or rounding source definitions.

For `SURFACE_OF_REVOLUTION`, the reader follows the STEP swept-curve and `AXIS1_PLACEMENT` fields.
An omitted axis direction defaults to z.
Four rational quadratic spans rotate each meridian control through one turn.
The rotation weights multiply the meridian weights, and the meridian knots remain unchanged.
A positive rational rescaling sets the largest axis component magnitude to one.
The four angular spans then stay within 120 degrees, and each rotation weight stays between one and four.
This avoids rounding an irrational unit axis or introducing extreme weights from an arbitrary direction magnitude.
The angular coordinate is a positive reparameterization of a full turn, so the source surface orientation is preserved.
The generating curve's parameterization remains the surface's v coordinate.
Control generation, weight multiplication and knot copying use the caller's allowance and deadline.
This exact carrier conversion supplies no face or solid certificate.

Solid construction separately certifies supported finite edges, complete embedded patch boundaries, topology and self-intersection.
General spline trims and unproved supports remain unsupported certification coverage.

### Measured source incidence

`measure_incidence` measures the entities named by a saved first-incidence refusal against the unchanged source.
Endpoint checks enclose the largest endpoint distance to the retained source curve.
Analytic support checks enclose the largest curve-to-support distance over the finite source feature when its complete parameter interval is known, or over an explicitly labelled enclosing carrier otherwise.
These quantities do not certify distance to a trimmed face or exported output.
Bernstein control hulls bound whole parameter cells; evaluated points supply lower witnesses.
The result records the attained bounds, their scope and whether the measurement width converged.

After an exact incidence refusal in a spline-bearing solid, `build_solid` may construct an exact planar derivative. The file entry also checks analytic sibling solids independently when the file contains splines.
This path requires planar supports with certified straight boundaries after repair. An exact finite-segment proof checks simple outer boundaries, disjoint holes, orientation and positive feature separation.
Shared vertices are projected to the exact common intersection of their adjacent planes.
The repair refuses new vertex mergers or unproved identity between distinct source edges.
Each shared spline is projected to its adjacent planes with unchanged knots and positive weights.
A positive line-image certificate permits a zero source parameter derivative at a span endpoint without admitting a reversal or constant span.
Exact homogeneous subdivision bounds its pointwise displacement against the complete retained source curve.
Line-segment and vertex displacement bounds are recorded separately.
Every attained bound must fit the existing 0.01 mm import cap and one tenth of the certified neighbouring feature size, including a conservative margin for the source-to-result displacement.
A separate certificate must bound each original vertex-to-curve gap within the import cap when both pieces move.
Planar projection/RREF, rational displacement-net subdivision, square-root enclosures, straight-image certificates and boundary predicates reserve the numeric layer's limb-based arithmetic cost before computation. Fraction comparisons reserve their cross products. These stages share the existing operation allowance and refuse when it is spent.
The derivative must pass the existing exact incidence, self-intersection and solid checks before admission.

Admitted bodies carry `import.incidence-repair` and `import.incidence-vertex-repair` distance marks for affected source entities.
Their JSON definitions name the solid and entity, the canonical source-geometry hash, attained deviation enclosure, approval cap and certified feature bound.
The source file and its recipe hash remain unchanged.
These local edge and vertex records are distinct from `face_tolerances` and do not claim a whole-face Hausdorff certificate.
Exact facts describe only the exact derivative.
An unsupported repair returns `incidence_gap` with its measured source bounds, missing certificate and any later strict solid refusal.

`catallax_mesh::preview_step` retains the exact source definitions and SHA-256 and returns an explicitly unproved `PreviewMesh` keyed by STEP entities, with shared edge samples and trimmed faces. It cannot supply a solid certificate, checked STEP or print mesh. Existing repair and recipe-publication gates are unchanged. See the [source and preview contract](../archive/interface-changes/2026-11-07-step-spline-source-preview.md) for supported encodings, refusal boundaries and verification.

### Tolerance-bearing source topology

`tolerant::read_step` retains the supported reader's source geometry in a shared vertex/edge graph without requiring exact incidence, loop order or solid closure.
The returned `TolerantFile` is source topology awaiting healing and validity checks; it grants no solid admission or exact-model extraction.
Each vertex and edge stores an individual numerical tolerance in millimetres.
`ToleranceOptions` defaults to an initial 0.0000001 mm and a maximum 0.01 mm; callers may state other nonnegative initial values under a positive cap before reading.
These controls describe local topological uncertainty, not a geometric movement or whole-body deviation certificate.

Growth follows [OCCT BRep_TVertex::UpdateTolerance](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingData/TKBRep/BRep/BRep_TVertex.lxx) and [BRep_TEdge::UpdateTolerance](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingData/TKBRep/BRep/BRep_TEdge.lxx): an entity's tolerance never decreases.
Every actual increase records its source entity, before/after values and named cause.
Growing an edge also grows either endpoint whose tolerance would otherwise be smaller, recording the triggering edge.
The complete update is prepaid against the active work allowance before its first write, and a cap or work failure preserves the previous graph.
An over-cap requirement refuses the update rather than clamping the recorded value below the requirement.
Any tolerance growth explicitly marks the body as repaired; zero growth establishes neither validity nor exactness.

The `catallax-import-tolerances/1` record saves controls, every current vertex/edge value and ordered growth history, scoped to the immutable source SHA-256 and body name.
`tolerant::reopen` rereads those exact bytes and reproduces the changes, including endpoint propagation, then compares the complete record.
Recognition retains the caller's Decider and deadline.
Reopen prepays the complete saved inventory and history, including quadratic work for variable-length numeric text, before parsing it.
Saved values use canonical integers or fractions; exponent notation is rejected before it can expand into a large integer.
The replay compares current entity values and history with explicit work charges, without constructing a second full record.
Changed or renumbered sources require separate correspondence and never inherit these bindings automatically.
This storage API does not add geometry-reading families or assign body tiers.
The existing certified and numerical import entry points retain their current behavior.

The default open shares one immutable Part 21 document between checked and native admission. Its header, SHA-256 and entity numbers stay tied to the source graph. Each route uses its existing unit rule when converting a recorded tolerance to millimetres.

### Transactional source healing

`tolerant::healing::heal` prepares a separate candidate from a `TolerantFile`.
Its stages are sewing, wire, face, shell and solid, in that order for each body.
The input graph and its immutable source SHA-256 remain available through `HealedFile::source`.
`definition` exposes the derived source definition; it does not construct an admitted solid.
`report` lists topology changes, measured local distances, checks actually passed and unresolved conditions.
Any topology change or tolerance growth sets `approximate`, including an orientation-only repair.
Later admission and STEP consumers must preserve that mark.

The stage order follows [OCCT ShapeFix_Shape](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingAlgorithms/TKShHealing/ShapeFix/ShapeFix_Shape.cxx), [ShapeFix_Wire](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingAlgorithms/TKShHealing/ShapeFix/ShapeFix_Wire.cxx) and [BRepBuilderAPI_Sewing](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingAlgorithms/TKTopAlgo/BRepBuilderAPI/BRepBuilderAPI_Sewing.cxx), with explicit candidate publication and the caller's cap.
The implemented repair family is conservative:

- Sewing pairs mutually unique free line edges on different faces of the same shell. Both finite endpoint correspondences must fit the cap; the final source-to-retained carrier projection checks also bound slightly skew line segments. Identical circle and ellipse carriers sew only with identical endpoints and the same finite traversal, including reversed endpoint order. Complementary arcs never match. A later vertex join that changes a sewn conic's finite domain refuses the complete candidate. Other nonlinear boundaries, ambiguous matches, different shells and different bodies remain separate.
- Wire connection first measures every declared adjacent gap and closure before applying endpoint joins. A disconnected ordering may instead use a unique walk through shared vertices, recording changed order and reversed uses. Small finite line intervals at or below the initial tolerance may be removed only when all incident supports are planar and every bound retains an edge. The measurement uses the original finite line, even after earlier joins. A closed circle is not a small line merely because its endpoints coincide.
- Planar line boundaries receive finite segment crossing/touch and orientation checks. Shell checks find connected components and consistent opposite edge uses, reporting open, nonmanifold or nonorientable cases. A closed planar candidate may have its whole orientation reversed from its signed polyhedral volume.

These operations retain source surfaces, line carriers, spline definitions and retained vertex coordinates.
They change explicitly recorded topology uses and orientation, rather than refitting source geometry.
Curved incidence, curved finite embedding, multiple-bound containment and cross-face validity remain named unresolved checks where the implementation cannot establish them.
Positive signed volume alone never grants solid validity.
A candidate with unresolved checks is available to the subsequent surface/solid tier consumer; this API itself assigns no tier.

Before publication, the healer checks every original vertex against its final retained target, every original line endpoint against its original carrier, and each original finite line interval against its incident source planes.
It also measures replaced finite line intervals against the retained carrier and removed small lines against their final retained vertex.
These final checks include source geometry removed from active uses, so sewing cannot hide a bad discarded carrier.
All required local tolerance maxima and edge-to-endpoint propagation are gathered and checked against the cap before ledger updates.
The ledger records committed increases; the separate measurements establish their geometric need.
These checks do not establish complete trimmed-face or whole-body source deviation.

Cap, topology, deadline or work refusal discards the complete private candidate, including earlier tentative repairs.
The caller's source graph and tolerance ledger remain unchanged.
`attempt` additionally returns the original typed refusal, failed stage/body and work spent on completed and partially failed stages, without exposing candidate geometry.
Preparation, all five stages and publication share the active allowance; a standalone call installs the unchanged 20M allowance.
Scans of cross-body change and issue logs prepay every potentially visited entry, including entries rejected by the body filter.
The final deadline checkpoint runs after the result fields are constructed.
`tolerant::reopen` also checks the deadline after its final saved-field comparison.

`HealedFile::tolerances` describes the derived graph's entity values and history.
The tolerance record alone cannot reconstruct topology changes; retaining the healing report and immutable source is required.
The report is an observation of these repairs, not an independently trusted solid certificate or a new persisted application format.
The existing recipe repair/approval path retains its own source-bound records and policy.
The [healer measurements](../archive/research/tolerant-import-healing/README.md) retain the actual fixtures and real-file outcomes.

### Per-body opening tiers and source replay

`tolerant::opened::open_step` first retains the immutable parsed graph, source units and body declarations. Each body is decoded only when its construction begins. `TolerantBody::declaration` identifies a deferred root; its decoded shell, vertex and edge tables are then empty, not a claim of empty geometry. `materialize` decodes that body under the caller's allowance without changing the retained declaration.
Each `BodyOutcome` reports `ValidSolid`, `RepairedSolid`, `TolerantSolid` or `SurfaceBody`, with an explicit approximation flag and actual refusal/check history.
`BodyView` distinguishes a checked model, a separately checked numerical B-rep, native sampled geometry and retained source surfaces.
Neither native sampled geometry nor retained source surfaces expose an exact `Model` or a conversion to certified solid authority.
An already admitted surface model retains its existing checked surface B-rep and surface kind.
The broad surface tier therefore does not itself assert finite embedding, closure or material; the view and recorded checks determine the available operations.

The default first runs the existing whole-file and strict per-body checked admission under the unchanged ordinary allowance. Successful bodies are settled and excluded from later native work. Refused or incomplete bodies then use the shared parsed source and cached carriers: main-compatible sampled native admission precedes the additional crossing, trace and pole repairs. Those additional repairs run only for bodies the checked paths could not admit.
Checked body retries decode the shared source through the existing analytic adapter before healing and numerical admission. Native fitting instead retains the original carrier graph. Deferring source decoding must not replace recognized checked-body geometry with unrecognized retained carriers.
After strict refusal, default source-import healing uses cap-certified point/line and line/plane incidence. The exact squared-distance comparison must prove correspondence within the caller's maximum. A nonzero incidence records that cap as a conservative proven upper bound, labelled `exact-squared-cap-upper-bound`; exact zero stays zero. It is not a measured residual or tightened deviation bound. Both candidate and original-source incidence checks retain their cap proofs. Sewing, wire repair, actual source-to-candidate movement checks, shell and solid checks, publication and native fitting are unchanged. Explicit checked opening, direct healing and parts built from scratch retain tightened bounds.
Both native passes validate an authored pcurve against the sampled edge and surface before retaining it. A first-pass mismatch refuses that candidate and sends the body to boundary repair; an unvalidated source trace cannot become an export fallback.
Wire ordering and connection, capped sewing, shell incidence, native curve/surface fits and boundary repairs operate on private body candidates. A refusal discards that candidate and its tentative evidence, not successful siblings. Native bodies retain stable `source_index` bindings across both passes. Optional stricter proof can upgrade a completed native candidate under the remaining ordinary allowance, but cannot discard it. A topology or 3D geometry repair remains approximate after an upgrade. Refused stitches retain the original source surfaces; a checked surface proof cannot publish the refused candidate as a solid.

After successful native fallback, `BodyOutcome::refusal` is empty. `checked_refusal` retains the earlier checked-tier failure separately, and `upgrade_refusal` records a later optional proof failure.
`NativeBodyFailure` records the first active repair/check stage, original typed refusal and body work. Bodies that cannot start after request exhaustion are marked `request_limit`. Resource refusals are not geometry diagnoses. `BodyOutcome::upgrade_refusal` is separate from a native-solid blocker, and `OpenedStep::tolerant_admission` exposes the native ledger even after partial admission or successful higher upgrades.

Missing or inconsistent pcurves are reconstructed from the 3D curve, with vertex incidence measured separately. Plane projections preserve B-spline knots and weights; analytic conics preserve their angular parameterization, including reflected UV frames. Other traces use sampled adaptive UV reconstruction at the initial precision and refuse unresolved reconstruction rather than claiming a complete deviation bound. Wire reconnection first tries a unique existing vertex walk. Adjacent finite-curve crossings can be cut at a uniquely measured intersection within the sewing/repair cap; ambiguous crossings refuse. This is not a general nonadjacent-wire arrangement algorithm. Missing pole edges are inserted only after measuring a collapsed surface-parameter gap, with explicit topology, pcurve and repair records.

Crossing candidates must identify the same finite parameter pair; equal XYZ at different parameters refuses rather than choosing a retraced subpath. Pole alignment remeasures changed final UV endpoints against both the 3D curve and topology vertex. Samples retain the final curve residual and a separate endpoint vertex residual; entity tolerances include both incidences. Recognized line previews use repaired endpoints in millimetres, not parameters of the original scaled STEP vector.

Sphere charts bound latitude by a 96-bit interval enclosure of the angular quarter-period, matching the [canonical spherical parameter domain](https://occt3d.com/dev/doc/refman/html/class_geom___spherical_surface.html). Latitude is not a periodic sheet. At a collapsed chart direction, native projection can step along the remaining live parameter direction; the existing residual-decreasing line search and convergence threshold still decide acceptance. This prevents a bad source trace from rebuilding on an unrelated latitude branch or leaving projection stuck at a pole.

Finite self-loop curves must have measured endpoint incidence and close within the initial entity tolerance before admission or sewing; a repeated topological vertex is not closure evidence. All assigned 3D angular contexts are resolved before interpreting source geometry. Disagreeing radians/degrees contexts refuse with `import.tolerant-ambiguous-units`; per-representation angular evaluation and export are not implemented.

`TolerantSolid` means these recorded local and topology checks passed. It does not establish whole-domain embedding, global self-intersection absence, outward material orientation or a complete source-to-body deviation bound. Natural edgeless faces must close both parameter directions, including collapsed pole rims; one closed spline direction does not cap an open tube. `open_checked_step` retains the explicit former certified/numerical policy.

`TolerantImport` exposes immutable source declarations, native curve/surface evaluators, sampled UV trims, entity tolerances, repairs and checks. `catallax_mesh::preview_tolerant_body` samples those carriers into an explicitly unproved display mesh. It samples through `DisplaySurface` and `DisplayCurve`, which evaluate the carrier's charts in binary64 with a running rounding bound and enclose only a sample binary64 cannot place within 1e-8 mm; admission and fitting keep the enclosing `NativeSurface` and `NativeCurve`. On a face with a retained native surface its chord and angle targets are measured on the carrier, so a vertex admitted off a face is drawn where it was admitted and does not drive refinement; a boundary sample farther from the carrier than the chord target plus the tolerance recorded for its edge and end vertices refuses the preview. Native STEP export rebuilds topology and preserves source carrier definitions, labels, schemas and the source occurrence graph. Whole-document export refuses incomplete native admission. Selected-body export returns a `NativeStepFile` containing body-local bytes, source-root bindings and export checks, with a fresh product frame rather than assembly placement. Shared source edges receive separate emitted topology and pcurve associations for each body; wrapped 3D bases cannot pull in sibling-only surface associations. Consistent nonlinear source pcurves are retained; reconstructed traces carry their measured correspondence and parameter conversion.

Native files have a short approximate `FILE_DESCRIPTION`. Concatenate the payloads of the ordered `CATALLAX_REPORT_JSON` comments and decode JSON to recover the repair/check record; JSON slash escapes prevent source strings from terminating those ignored STEP comments. Body-local exports include only the selected body's geometry, tolerances, repairs and checks, plus source-wide authority limitations.

`OpenedStep::export_step_document` writes the default selector's actual certified, numerical or native body views, without another admission. The selected body artifacts keep their own dimensional contexts and are referenced through representation maps from the original prototype bindings. Source occurrence links and placement transforms remain attached to those bindings. The result records original prototype identities, emitted geometry roots and per-body authority separately; it does not promote an approximate body to certified authority.

Whole-source selected export requires a retained parsed graph and every source prototype reachable from a product-bound shape representation. A free representation does not establish product or occurrence coverage. Any unreachable prototype or surface-only sibling refuses the document atomically. Body writers, source-graph composition and serialization use the caller's existing export allowance. The public `repeated_placements` control covers one shared part at three distinct transforms in millimetre and metre source contexts, plus a valid but detached prototype that must refuse whole-document export.

Normal default reopening resolves prototype geometry units from the representations that directly contain those prototypes. Assembly placement contexts can use different units without reinterpreting body coordinates. All bound prototype contexts must still agree on length and plane-angle units; conflicting contexts for the same prototype refuse. Any unbound prototype retains the globally unambiguous-unit requirement. Native body-local export selects a context matching its compiled geometry units rather than the first source context. This does not add admission of prototypes with different geometry units.

Default controls remain 0.0000001 mm initial entity tolerance and a 0.01 mm sewing displacement cap. Native measured entity tolerances may grow beyond that cap and are reported individually; this is not the fixed-cap checked policy. The library introduces no implicit native work allowance. The existing census observer still applies a shared 20M opening allowance and records exhausted work separately from geometry refusals.
Output slots are reserved before optional admission work.
A body-local budget/deadline refusal discards that private candidate and leaves its immutable source declaration. Completed siblings and native bodies awaiting optional proof remain usable.
An `Ok` result can therefore contain only retained surfaces even when optional work exhausted the outer allowance.
Callers must handle body views explicitly instead of interpreting success as solid admission.
Source hash mismatch, unreadable document/schema data and ambiguous geometry units remain document-level typed failures. A readable body declaration can remain retained when its own topology or geometry construction refuses.
Existing recipe construction and proposal/approval entries retain their explicit checked contracts; this default changes the opened-import consumer, not certified model construction.

`OpenedStep::saved_controls` saves the original `catallax-import-tolerances/1` source record.
`opened::reopen` verifies that record against the immutable source, rejects a purported source record containing later tolerance growth, then reruns the selector with the recorded controls and the caller's step identity.
It introduces no healed-body codec or cached-model promotion.
The replay and reconstruction share one active allowance.
The normal 0.001 mm gap fixture saves/reloads its controls and reproduces an approximate six-face, twelve-edge, eight-vertex body with identical repair history.
All entity keys remain scoped to the original source SHA-256 and body; changing or renumbering a file still requires the separate verified correspondence/replacement path.

The census observer separates supported source retention, admitted solid tiers and independently checked exports.
Its source inventory coverage flag is the frozen source-inspection result at the identical source hash; incomplete representation coverage is never promoted to a complete-source solid file.
Retained source surfaces without an admitted B-rep receive no checked STEP or valid-solid credit.
All admitted bodies share one separate 20M export allowance, and actual written bytes are independently checked with the pinned OCCT builder protocol.
The tier collector also caps each producer and independent reader at 6,144 MiB by default, using the shared RSS guard. Timeouts, memory terminations and process failures remain unmeasured rows rather than geometric losses. An unmeasured baseline pass prevents acceptance.
A producer crash after its complete opening checkpoint also retains that opening verdict; only the unfinished export is a process error.
The observer separately reports files with at least one individually valid exported solid and strict complete-source/all-export coverage. Strict source coverage requires an occurrence-preserving whole-document export, including when every prototype obtains a higher checked tier. Isolated prototype exports do not establish source occurrence coverage. Historical prototype-valid totals remain a separate comparison. A valid sibling counts for existence even if another root fails. Incomplete export observation makes the observed count a lower bound; zero observed solids is then unknown, not a negative verdict on skipped siblings. Raw-root validity, emitted root topology counts and default occurrence-transfer validity are separate observations. Neither these checks nor the OCCT reference's default transfer plus ShapeFix establishes a common tolerance or source-fidelity result.
The independent reader also traverses the emitted product-bound representation graph and checks every claimed body root against it. A valid detached body can receive prototype/existence credit, but not strict complete-source credit; a generic successful root transfer does not override missing occurrence reachability.
Builder processing can heal internally; neither unhealed validity nor complete source fidelity follows from that result.
Each native root records expected and inspected parts. A converter exception can leave later parts uninspected, so completed producer export does not establish completed raw validation. Such roots, raw-reader timeout/OOM and missing results appear as unknown or observed lower bounds. Non-budget whole-document export refusals before the body loop are counted separately before considering any selected-body fallback.
The [frozen tier census](../archive/research/step-census-opened/README.md) owns the measured counts, producer and limitations.

`collect --tier-baseline PRIOR.csv` additionally pins a hash-matched prior census for per-file tier comparison. Recorded tier totals are compared separately from solid coverage, including co-occurring checked-count loss and tolerant-count gain while a file stays full. Missing tier inventories or changed body inventories produce explicit unknown comparisons. Historical aggregate totals establish neither body-identity correspondence nor a lower bound on same-body demotions.

### IGES source recognition

`read_iges` and `read_iges_with_decider` read IGES 5.3 manifold solid B-rep objects into the same `StepFile`, so `build_file`, its exact incidence checks, the self-intersection check and every builder refusal are shared. `import_iges_bodies` checks the saved hash, reads the file, checks self-intersection and builds, as `import_step_bodies` does for STEP, including its separately recorded bounded spline-incidence path. The recipe step `import_iges` (RECIPE.md section 4.18.2) goes through `apply_source` with `SourceFormat::Iges`, so gap proposals (`propose_source`), saved decisions, repairs, recorded face tolerances and the STEP derivative are those of a STEP import. A persisted IGES record replays through `rebuild_source` with `SourceFormat::Iges`; `rebuild` is its STEP form, so no replay guesses the format from the bytes or their name. Face keys are directory-entry sequence numbers such as `D17`; edges and vertices are named by their list and index, such as `D3:5`.

Read: 186 with no voids, 514 form 1, 510, 508, 504 and 502 form 1; surfaces 190, 192, 194, 196 and 198 in forms 0 and 1, and 128; curves 110, 100 and 126; 116 points, 123 directions and 124 transforms, chained up to eight. Loops run with the face on their left about the surface normal, so a face against its surface (shell flag 0, or a 186 flag 0) reverses them. Parameter-space curves in loops are not read; the model-space curve decides. A vertex use among a loop's edges is a degenerate edge, such as a sphere pole between two uses of its seam, and bounds nothing. A two-point degree-one 126 with equal weights is a line. A 194 at exactly 45 degrees is the exact 45-degree cone; any other angle is a degree semi-angle, built through its rims as for STEP. Lengths convert to millimetres exactly from the unit flag or name (inch, millimetre, foot, mile, metre, kilometre, mil, micron, centimetre, microinch); a model-space scale other than 1 is `UnknownUnit`.

Typed refusals: `Unsupported` names the entity and IGES type for other entities and forms, a 186 with voids, an open shell, a transform that is not exactly orthonormal, a topology entity with a transform, and a 100 arc whose radius is not rational. A file with no 186 names its first 143, 144 or 514 entity as `Unsupported`, or is `Malformed`. Non-ASCII text is `Malformed` before any column is read, and every count in the file is checked against the parameters its entity holds, with checked arithmetic, before anything is allocated or indexed from it. Files from other writers often carry binary trigonometric noise: an OpenCASCADE cylinder's seam vertex at `y = -7.3e-16` is `not_solid`, and a matrix entry of `2.4e-16` is `Unsupported`. Neither is moved; curved repairs with recorded local tolerances are separate work. See [the interface note](../archive/interface-changes/2026-11-11-iges-import.md).

### Binary-noise cleanup

Files written from binary floating point put a seam vertex at `y = -7.3e-16` or a rotation entry at `2.4e-16` where the writer meant zero, so their exact reading refuses: a vertex off its curve (`not_solid`), a rotation that is not exactly orthonormal (`Unsupported`), or a self-intersection. `apply_source` (and so `apply`, `rebuild_source` and the recipe import steps) then tries a cleaned reading of the same bytes. Every decimal written with at least ten significant digits is replaced by the simplest rational within a relative `1e-12` of it when that rational's denominator is at most 10,000; shorter decimals are taken as written, and values that are not near a simple rational, such as `0.7071067811865476`, are unchanged. The cleaned file must then pass every exact check of an ordinary import: incidence, self-intersection, gap policy and the scratch rebuild.

The movement is proven for the whole file. A changed length moves a face by at most its change; a changed ratio turns a direction whose squared norm is at least a quarter of the file's shortest by at most eight times the change over that quarter, across at most four times the largest magnitude read. The sum must be at most the 0.01 mm automatic cap and a tenth of the file's conservative feature size (the shortest edge chord or radius). Files with spline data are not cleaned. Every face is then listed in `face_tolerances` with that outward-rounded upper bound, a lower bound of zero and the automatic approval, and the cleaned file is the derivative, which `rebuild_source` reproduces. Otherwise the exact refusal stands. The exact reading always comes first, so a file that imports exactly is never cleaned.

## Manufacturing information

`read_step_pmi(bytes)` reads a file's AP242 semantic PMI into a
`catallax_step::pmi::Pmi<String>`. Faces are named by the entity of their
`ADVANCED_FACE`, such as `#17`: the face keys an import of the same bytes
reports. It reads what [piece 7 writes](07-step-export.md#manufacturing-information)
and the common forms of the CAx-IF recommended practices, in any length unit.

A dimension's nominal is the value the file states, in millimetres. The reader
does not compare it with the geometry; the caller measures the imported body
with `catallax_step::pmi::nominals` and decides.

A dimension with the `theoretical` note reads as basic (`Limits::Basic`). A
diameter on a feature that carries the attribute `hole depth` with the value
`through` reads with `Size::through` set.

Nothing is guessed. An item the reader cannot represent in full is left out
and listed in `StepPmi::skipped` with its entity, type and reason: an angular
dimension, a dimension along a path, a reference dimension or one with
another modifier or note, a datum target other than a placed
target on one planar face, a common datum, a tolerance with a modifier other
than a material condition, a datum reference with a modifier other than a
material boundary, an item on geometry other than faces. A datum that only
datum targets set up, with no `DATUM_FEATURE` naming its faces, is skipped
with its targets. A composite frame's lower segment is skipped with its upper. A
tolerance is skipped if no accepted datum defines the letter it refers to.
A datum skipped because its faces already carry a datum does not reserve
its letter. A later datum on a different feature can define that letter.
The only refusal is a file that is not readable Part 21 text.

`cargo run -p catallax-import --example step_pmi -- part.step` prints what a
file reads as and every item left out.

## Immutable source

The file bytes are stored unchanged and addressed by `sha256`. The recipe step keeps `file` and `decisions`. A repaired derivative is a separate artifact derived from (source hash, approved change record); it is never written back over the source and can be rebuilt from the two inputs. A hash mismatch is a recipe error, not geometry.

## Repair policy

Measurement comes first. A gap is a pair of boundary features (edge/edge or vertex/vertex) that the file declares as distinct but that are within distance `d` of each other. Distance is computed exactly or with a certified enclosure; a bound that cannot decide against a threshold counts as uncertain.

| Case | Action |
|---|---|
| `d <= 0.01 mm`, both features on one shell, closing joins two faces of that shell, and the move is below one tenth of the nearest-feature size | Repair automatically, record it |
| Features on different shells or bodies, any `d` | Never merge. Bodies stay separate, even at 0.005 mm |
| The gap is an opening (a boundary loop of a face with no neighbouring face, an open shell, a hole) | Never close automatically, even at 0.005 mm. Ask |
| `d > 0.01 mm`, or `d` or the nearest-feature size is uncertain, or identity or opening intent cannot be established | Ask. Proximity and any AI suggestion are not permission |
| Declined | The source bytes are never modified and no derivative is made. The decision is saved as `decline` so the question is not re-asked for the same bytes. The solid that has the declined gap stays open and is not built; it is reported in `left_open`, and the application must tell the user so. This is a documented limit, not a returned original body |

An approved out-of-policy proposal is still validated deterministically after the move: closed, non-self-intersecting, oriented, connected, and the measured maximum movement must not exceed the approved tolerance. Failed validation rejects the repair and leaves the original.

## Change record (`catallax-import-repair/2`)

The recipe `decisions` saves answers under the proposal's exact `question`, `repair:<sha256-prefix>:g<gap-index>`. The JSON change record is kept beside the recipe so reopening restores it.

```
{
  "format": "catallax-import-repair/2",
  "source": {"name": "...", "sha256": "..."},
  "derivative": {"sha256": "..."},
  "changes": [{
    "id": "g0",
    "entities": ["#12", "#40"],
    "kind": "seam_close",
    "bounds_before": [[x0,y0,z0],[x1,y1,z1]],
    "bounds_after":  [[x0,y0,z0],[x1,y1,z1]],
    "max_movement_mm": "0.004",
    "tolerance_mm": "0.01",
    "reason": "boundary gap below automatic cap",
    "approval": "automatic" | "user" | "declined"
  }],
  "face_tolerances": [{
    "solid": "#178",
    "face": "#40",
    "deviation_lower_mm": "0.004",
    "deviation_upper_mm": "0.004",
    "approved_tolerance_mm": "0.01",
    "feature_size_squared_mm2": "100"
  }]
}
```

Distance numbers are decimal strings; `feature_size_squared_mm2` is an exact rational string in square millimetres. `max_movement_mm` measures repaired vertex displacement, not whole-face deviation. A `user` approval records the approver identity the app supplies and the approved tolerance. `face_tolerances` lists only nonzero changed-face deviations and rim-built cones; it is empty for clean imports without such cones and for declined solids. Rebuild verifies every bound and the derivative hash. Old `/1` records are not accepted as `/2` evidence: reapply the saved source and decisions, then persist the newly verified record.

## Acceptance tests (public synthetic fixtures)

1. A box whose two faces miss by 0.004 mm at one seam repairs; movement at most 0.004 mm; record written.
2. Two boxes 0.005 mm apart import as two bodies; no record of merging.
3. A box with a 0.005 mm slot through one wall stays open; the question is asked.
4. A 0.02 mm gap and an undecidable gap each ask; decline returns bytes-identical original geometry and no derivative.
5. Approve an over-cap gap: record survives save and reopen, and rebuild from (source, record) reproduces the derivative hash.
6. Unsupported entity: refusal names the entity; nothing is partially imported.
7. Controls: 1 and 2 and 4 workers, cold and warm cache, a scaled and transformed copy, and exact ties at 0.01 mm (accepted) and just above it (asked). Existing 1 nm, 100 nm, 1 micrometre and tenth-of-nearest-feature limits stay unchanged.

Tier when built: FULL SIX for interpretation and numerics; dialog-only changes SINGLEPAIR.

## File ownership needed

`catallax-recipe` (consumer of `ImportStep`), `catallax-topology` (shell construction), `catallax-app` adapter entry points (L03), cache mechanics (L04), and a new `catallax-step` reader module. SolidDraft: import/repair dialog, change-record viewer, decline path.

## Implemented repair scope

`propose` and `apply` in `catallax-import` repair straight planar incidence as described above and measure straight seams between planar faces as follows.

- A gap is a set of straight boundary edges used once with no opposite use, paired within a search radius (default 1 mm) by opposite travel and endpoint distance. Pairs sharing a point form one gap. Distances are exact squared rational distances; the shown bounds enclose their square roots at 1e-12 mm.
- A gap is an opening, always asked and never closed, when any of its boundary edges has no partner within the radius, or when a pair's faces are exactly anti-parallel (closing would fold them onto each other). A thin plate with one edge face missing is such an opening. Accepting an opening proposal is refused; only a decline is accepted.
- Automatic closure needs a seam, a distance of at most 0.01 mm inclusive and a squared distance times 100 not above the conservative squared feature size of the faces involved. Feature size includes straight-edge lengths and nonincident vertex distances to supporting edge lines of proven convex polygons. Missing feature evidence asks. Final whole-face bounds must also meet the feature limit, even with user approval.
- The repair moves whole vertices of one side onto their partners and refits each planar face through its moved vertices, which must stay exactly coplanar. Both sides are tried. The side that changes fewer surfaces wins, then the smaller movement. So a seam closes by extending faces along their own surfaces where possible. If neither side works, `validation_failed` and nothing is built.
- The repaired solid is rebuilt in a scratch model, checked for closed, sound topology and positive volume before it is built for the caller. Self-intersection is tested exactly (see below) on the exact-incidence source or the certified incidence candidate and again on every repaired solid.
- Recorded body bounds cover the solid's vertices; they ignore curved bulge. Movement is the largest vertex displacement, rounded up at 1e-12 mm. Separate source-relative whole-face distance records use the convex-cover proof above and the same outward rounding.
- Curved gaps are out of scope for the first path. A gap that involves a curved edge or face is refused with `not_solid` or `validation_failed`, never guessed at. The message names the faces and edges involved, and the gap id once one is assigned. This is a documented limit.
- A decline leaves that solid open and unbuilt: it is reported in `left_open` and not built, the source bytes are never modified, and the application must say plainly that the solid was left open. This is a documented limit.
- A solid with no gap or incidence repair is built as written. Separate solids are never joined, so bodies 0.005 mm apart produce no proposal.

### Self-intersection

`check_self_intersection` and the build entries compare every pair of faces of one solid and refuse with `self_intersecting`, naming the solid and the first pair of faces, or with `self_intersection_not_verified` when a pair cannot be decided. A pass means every pair was decided. Separate solids are never compared with each other.

The rule is the same for every kind of face. Two faces intersect when their carrier surfaces meet in a curve and a positive-length piece of it lies in both closed faces and in the interior of at least one; or when the carriers coincide, the faces look the same way and their interiors overlap. Contact along a shared edge, at a single point (including a sphere tangent to a face) and back-to-back coincident faces facing opposite ways are allowed. A cylinder wall that meets its own cap along their common circle passes. A boss whose ring is a hole boundary of the face it stands on passes. A boss standing on a whole face with no hole for it is refused, because its ring lies in that face's interior.

Planar faces with straight edges are compared exactly in rational arithmetic. Every other pair (cylinder, 45 degree cone, sphere, torus and any face with a circular or elliptical edge) goes to the certified intersection engine, after pairs whose proven boxes are strictly apart are set aside. The engine returns the curves where the carriers meet. Each curve is cut where it meets the edges of either face, one exact point inside every piece is found, and its membership in both faces and on their boundaries is decided with exact or certified answers. No floating-point value decides any step. An engine answer that is not exact is treated as undecided.

Not verified, with the pair and the reason named, and never a pass: a face with no boundary on a plane, cylinder or cone; coincident carriers facing the same way where a region is bounded by curved edges; a boundary edge that lies in both carriers and meets the meeting curve without matching it; any pair the engine refuses or cannot make exact; and a solid whose faces cannot be set up for the engine. An intersection found for any pair wins over an undecided pair. A solid with an undecided pair is refused, so it cannot be imported until the pair is decidable.

Not fixtured: a repair that itself creates an intersection. The repaired path shares the checker, and no fixture forces such a repair.

## Selecting one source solid

`ImportStep` and `ImportIges` accept an optional `source_solid` string. Call
`catallax_import::SourceFormat::solid_entities(bytes)` to list the keys and use
its returned length as the solid count. STEP keys are solid entity references
such as `#150`; IGES keys are solid directory references such as `D301`.
Store a key with the original file's SHA-256. It identifies that source entity,
including when several solids share faces. It is neither a face-set fingerprint
nor an ordinal, and is not an identity across changed source bytes.

```json
{
  "id": "body_one",
  "op": "import_step",
  "file": {"name": "reviewed.step", "sha256": "<original source SHA-256>"},
  "source_solid": "#150",
  "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}
}
```

The evaluator runs the existing source import and repair decisions, then selects
that entity before collecting solids, placing the result and creating its
checkpoint. `Built::import_solids(step)` contains only the selected entity and
its placed lumps. Face origins keep their original source keys. Repair evidence
still covers the original file, including decisions about its other solids.
An unrelated declined solid does not prevent selection of a closed solid, but
missing required repair decisions and invalid source geometry still refuse
through the shared importer. Selecting a declined solid returns
`import.left-open`; an unknown key returns `import.source-solid`.

Applications can build one recipe and `EvaluationCheckpoint` per listed key
using the same original bytes, hash and reviewed decisions. Pass each checkpoint
to `evaluate_incremental_with_files` when adding downstream features. Changing
the selector invalidates prefix reuse. Omitting the selector retains the existing
whole-file behavior and omits the field from serialized recipe bytes. Rust callers
constructing either operation directly must supply `source_solid: None` to retain
that behavior.

`import_edit::replace_step` refuses a hash-changing replacement of a selected
import with `import.replace-source-solid`. Its replacement proposal has no
explicit new solid selector, so carrying an old entity key into different bytes
could select another solid. Rebuilding with the same source hash remains allowed.
