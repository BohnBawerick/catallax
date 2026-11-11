# Piece 8: STEP import and bounded gap repair (contract)

Status: `catallax-import` implements the reader, body construction, repair policy and change records described below. The file-carrying feature and application build entries evaluate `Operation::ImportStep` for one solid, with saved approvals and declines, retained source-relative face tolerances and immutable source verification. SolidDraft preview, persistence and checked-export integration remain separate companion work. See `docs/interface-changes/2026-11-02-planar-import-local-tolerances.md` for the current build and record contract.

## Application roundtrip slice

`crates/catallax-app/tests/import_roundtrip.rs` exports a 2 × 3 × 5 mm block and a radius-3, height-7 mm cylinder with `catallax-step::write_step`, then calls `catallax-import::propose`, `apply` and `rebuild` and `catallax_app::build_incremental_with_files`. These entries already existed; this slice adds application-level acceptance evidence, not another importer or a new wire format.

Both sources need no repairs and produce no derivative. Exact carrier, edge-curve and vertex coincidence queries compare the application body with the source primitives. The checks also require sound topology, the file's face origins, matching exact bounding boxes, overlapping volume enclosures narrower than 0.000001 mm³, and a valid checked STEP re-export. Import movement is zero for these fixtures. The source export's reported bound remains a separate quantity and is not inferred from volume agreement.

Display is checked through the application entry at kernel defaults and at SolidDraft's automatic coarse, standard and fine controls, with chord height `0.002 × body diagonal × {4, 1, 0.25}` mm and a 15-degree surface-normal bound. Every face must meet both requested bounds within the existing application work guard. This tests Catallax's controls, not SolidDraft's f32 conversion, rendered viewport or complete application workflow. No outside STEP fixtures are used.

Still refused or unproven by the clean roundtrip slice: splines, curved gap repair and multiple solids in one recipe import step. Spheres, tori and cones of any semi-angle now build (see the surface list below), but this slice does not exercise them. The lower-level importer can return separate bodies; the recipe application does not merge them. Seam exchange contracts and field construction from imported faces are unchanged.

SolidDraft uses `propose` for the preview and `apply` for the repair record, saves unchanged source bytes and their SHA-256 plus the recipe decisions, and supplies those bytes through `ImportFiles` to `build_incremental_with_files`. It must use each proposal's exact `question` string for saved decisions, not construct its own gap-index spelling. `rebuild` checks the persisted repair record on reopen. Recipe values are `decline` or `accept:<tolerance-mm>:<approver>`. Automatic repairs save no acceptance; user acceptance applies only to asked seam proposals.

## Nonzero planar repair qualification

`crates/catallax-app/tests/import_step.rs` exercises nonzero repairs through proposal, application, serialized recipe and repair record, rebuild, cold application rebuild, retained-checkpoint display and checked STEP output. Cases include a 0.004 mm automatic repair, the inclusive 0.01 mm cap, a just-over-cap saved approval, a 0.02 mm saved approval, saved declines, and the inclusive tenth-of-feature limit with a just-over-limit approval refused.

Repair records now use `catallax-import-repair/2`. Each changed convex planar face has an enclosure of its two-sided Hausdorff distance to the complete source trimmed face. The proof checks straight-edge incidence, planarity, ordered strict convexity and complete polygon coverage. Exact vertex-to-trim distances provide lower witnesses. Distance to a convex polygon is convex, so the same witnesses bound every triangle in a complete fan in both directions. This is not a relabelled vertex displacement or carrier-plane residual. Faces with holes, nonconvex covers or unproved adjacent feature bounds refuse rather than acquire a certificate.

The attained upper bound is the local tolerance. It must fit the approved file limit and a tenth of a conservative nearest-feature bound, including narrow skew faces and source/repaired neighbours. A user can approve a larger file limit, not override the feature limit. Bounds are remeasured against the original source after all repairs, never summed or clamped.

`Built::import_record(step)` retains the source-bound record with its immutable import body and recipe, including checkpoint reuse. Rigid placement preserves distances. These records do not certify later modified faces, a whole solid against an intended design, display rounding or exported output. STEP's existing error report remains relative to the repaired model. To bound output against a source face, a consumer needs the applicable source-face record, output correspondence and the separate output bound.

The standalone persisted-artifact smoke approved a 0.02 mm seam and reread the saved recipe, unchanged source, `/2` repair record and derivative. It produced four exact `[0.02, 0.02]` mm face enclosures, a six-face, twelve-triangle OBJ and a checked STEP re-export with zero output error against the repaired model. The separate approval regression grants 0.03 mm while requiring the attained face bounds to remain 0.02 mm. These are API and artifact checks, not SolidDraft viewport or independent OpenCASCADE validation.

## Supported first path

Files: AP203/AP214/AP242 `MANIFOLD_SOLID_BREP` with one `CLOSED_SHELL`, in millimetres, inches or metres (exact conversion).

Surfaces: `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`. Curves: `LINE`, `CIRCLE`, `ELLIPSE`, and `SURFACE_CURVE` or `SEAM_CURVE` read through their 3D curve (parameter-space curves are ignored). A `CONICAL_SURFACE` at exactly 45 degrees in a degree context is built exactly as declared. At any other semi-angle, in degrees or radians, tan(semi-angle) is irrational, so no decimal rim away from the placement plane lies exactly on the declared cone. Such a cone is built through the face's two coaxial rim circles, or through its rim and its apex. The apex is a `VERTEX_LOOP`, or the end of a seam on the axis. That cone has a rational slope. Its whole-face deviation from the declared cone is certified. The largest radial miss at the two ends gives the upper bound. The lower bound is an end's distance to the declared double cone: the smaller of its misses from the two nappes' lines in the meridian plane, times cos(semi-angle). The cone is accepted automatically only when the upper bound fits the 0.01 mm automatic cap and a tenth of the face's nearest feature. The feature is the axial span or a nonzero rim radius. The deviation is recorded in `face_tolerances`. Beyond either limit the file is refused as `not_solid`, with both bounds in the detail. A sphere or torus face whose bounds reduce to seams, or to one `VERTEX_LOOP` on the surface, is the complete surface. Circles on spheres and torus parallels and meridians are checked exactly. A `VERTEX_LOOP` elsewhere is refused. Several `MANIFOLD_SOLID_BREP` entities import as separate bodies, which acceptance test 2 needs. Offset, swept, voids, shell models and assemblies are refused with the entity id and type. A face with a single `FACE_BOUND` and no `FACE_OUTER_BOUND` is read with that bound as its outer boundary, as Catallax's own exporter writes it. Decimals read as exact rationals; no float enters a decision.

Face keys are the `ADVANCED_FACE` entity names, as in `docs/RECIPE.md` section 4.18.

### Spline source recognition and preview

`read_step` and `read_step_with_decider` additionally read simple and complex knot-bearing B-spline curves and surfaces, positive rational weights, and Bezier, uniform and quasi-uniform encodings. Exact rational knot extraction clamps the active domain without fitting or rounding source definitions. This is source recognition, not certified solid construction: `build_file` and `build_solid` refuse general spline faces and edges as unsupported certification coverage.

`catallax_mesh::preview_step` retains the exact source definitions and SHA-256 and returns an explicitly unproved `PreviewMesh` keyed by STEP entities, with shared edge samples and trimmed faces. It cannot supply a solid certificate, checked STEP or print mesh. Existing repair and recipe-publication gates are unchanged. See the [source and preview contract](../interface-changes/2026-11-07-step-spline-source-preview.md) for supported encodings, refusal boundaries and verification.

### IGES source recognition

`read_iges` and `read_iges_with_decider` read IGES 5.3 manifold solid B-rep objects into the same `StepFile`, so `build_file`, its exact incidence checks, the self-intersection check and every builder refusal are shared. `import_iges_bodies` checks the saved hash, reads the file, checks self-intersection and builds, as `import_step_bodies` does for STEP; nothing is repaired. The recipe step `import_iges` (RECIPE.md section 4.18.2) goes through `apply_source` with `SourceFormat::Iges`, so gap proposals (`propose_source`), saved decisions, repairs, recorded face tolerances and the STEP derivative are those of a STEP import. A persisted IGES record replays through `rebuild_source` with `SourceFormat::Iges`; `rebuild` is its STEP form, so no replay guesses the format from the bytes or their name. Face keys are directory-entry sequence numbers such as `D17`; edges and vertices are named by their list and index, such as `D3:5`.

Read: 186 with no voids, 514 form 1, 510, 508, 504 and 502 form 1; surfaces 190, 192, 194, 196 and 198 in forms 0 and 1, and 128; curves 110, 100 and 126; 116 points, 123 directions and 124 transforms, chained up to eight. Loops run with the face on their left about the surface normal, so a face against its surface (shell flag 0, or a 186 flag 0) reverses them. Parameter-space curves in loops are not read; the model-space curve decides. A vertex use among a loop's edges is a degenerate edge, such as a sphere pole between two uses of its seam, and bounds nothing. A two-point degree-one 126 with equal weights is a line. A 194 at exactly 45 degrees is the exact 45-degree cone; any other angle is a degree semi-angle, built through its rims as for STEP. Lengths convert to millimetres exactly from the unit flag or name (inch, millimetre, foot, mile, metre, kilometre, mil, micron, centimetre, microinch); a model-space scale other than 1 is `UnknownUnit`.

Typed refusals: `Unsupported` names the entity and IGES type for other entities and forms, a 186 with voids, an open shell, a transform that is not exactly orthonormal, a topology entity with a transform, and a 100 arc whose radius is not rational. A file with no 186 names its first 143, 144 or 514 entity as `Unsupported`, or is `Malformed`. Non-ASCII text is `Malformed` before any column is read, and every count in the file is checked against the parameters its entity holds, with checked arithmetic, before anything is allocated or indexed from it. Files from other writers often carry binary trigonometric noise: an OpenCASCADE cylinder's seam vertex at `y = -7.3e-16` is `not_solid`, and a matrix entry of `2.4e-16` is `Unsupported`. Neither is moved; curved repairs with recorded local tolerances are separate work. See [the interface note](../interface-changes/2026-11-11-iges-import.md).

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

`propose` and `apply` in `catallax-import` measure and repair only straight seams between planar faces.

- A gap is a set of straight boundary edges used once with no opposite use, paired within a search radius (default 1 mm) by opposite travel and endpoint distance. Pairs sharing a point form one gap. Distances are exact squared rational distances; the shown bounds enclose their square roots at 1e-12 mm.
- A gap is an opening, always asked and never closed, when any of its boundary edges has no partner within the radius, or when a pair's faces are exactly anti-parallel (closing would fold them onto each other). A thin plate with one edge face missing is such an opening. Accepting an opening proposal is refused; only a decline is accepted.
- Automatic closure needs a seam, a distance of at most 0.01 mm inclusive and a squared distance times 100 not above the conservative squared feature size of the faces involved. Feature size includes straight-edge lengths and nonincident vertex distances to supporting edge lines of proven convex polygons. Missing feature evidence asks. Final whole-face bounds must also meet the feature limit, even with user approval.
- The repair moves whole vertices of one side onto their partners and refits each planar face through its moved vertices, which must stay exactly coplanar. Both sides are tried. The side that changes fewer surfaces wins, then the smaller movement. So a seam closes by extending faces along their own surfaces where possible. If neither side works, `validation_failed` and nothing is built.
- The repaired solid is rebuilt in a scratch model, checked for closed, sound topology and positive volume before it is built for the caller. Self-intersection is tested exactly (see below) on the source solids and again on every repaired solid.
- Recorded body bounds cover the solid's vertices; they ignore curved bulge. Movement is the largest vertex displacement, rounded up at 1e-12 mm. Separate source-relative whole-face distance records use the convex-cover proof above and the same outward rounding.
- Curved gaps are out of scope for the first path. A gap that involves a curved edge or face is refused with `not_solid` or `validation_failed`, never guessed at. The message names the faces and edges involved, and the gap id once one is assigned. This is a documented limit.
- A decline leaves that solid open and unbuilt: it is reported in `left_open` and not built, the source bytes are never modified, and the application must say plainly that the solid was left open. This is a documented limit.
- A solid with no gap is built as written. Separate solids are never joined, so bodies 0.005 mm apart produce no proposal.

### Self-intersection

`check_self_intersection` and the build entries compare every pair of faces of one solid and refuse with `self_intersecting`, naming the solid and the first pair of faces, or with `self_intersection_not_verified` when a pair cannot be decided. A pass means every pair was decided. Separate solids are never compared with each other.

The rule is the same for every kind of face. Two faces intersect when their carrier surfaces meet in a curve and a positive-length piece of it lies in both closed faces and in the interior of at least one; or when the carriers coincide, the faces look the same way and their interiors overlap. Contact along a shared edge, at a single point (including a sphere tangent to a face) and back-to-back coincident faces facing opposite ways are allowed. A cylinder wall that meets its own cap along their common circle passes. A boss whose ring is a hole boundary of the face it stands on passes. A boss standing on a whole face with no hole for it is refused, because its ring lies in that face's interior.

Planar faces with straight edges are compared exactly in rational arithmetic. Every other pair (cylinder, 45 degree cone, sphere, torus and any face with a circular or elliptical edge) goes to the certified intersection engine, after pairs whose proven boxes are strictly apart are set aside. The engine returns the curves where the carriers meet. Each curve is cut where it meets the edges of either face, one exact point inside every piece is found, and its membership in both faces and on their boundaries is decided with exact or certified answers. No floating-point value decides any step. An engine answer that is not exact is treated as undecided.

Not verified, with the pair and the reason named, and never a pass: a face with no boundary on a plane, cylinder or cone; coincident carriers facing the same way where a region is bounded by curved edges; a boundary edge that lies in both carriers and meets the meeting curve without matching it; any pair the engine refuses or cannot make exact; and a solid whose faces cannot be set up for the engine. An intersection found for any pair wins over an undecided pair. A solid with an undecided pair is refused, so it cannot be imported until the pair is decidable.

Not fixtured: a repair that itself creates an intersection. The repaired path shares the checker, and no fixture forces such a repair.
