# Piece 8: STEP import and bounded gap repair (contract)

Status: `catallax-import` implements the reader, body construction, repair policy and change records described below. The file-carrying feature and application build entries evaluate `Operation::ImportStep` for one solid, with saved declines and immutable source verification. SolidDraft preview, persistence and checked-export integration remain the separate V9-S companion. See `docs/interface-changes/2026-10-30-step-import-recipe-evaluation.md` for the build contract.

## Application roundtrip slice

`crates/catallax-app/tests/import_roundtrip.rs` exports a 2 × 3 × 5 mm block and a radius-3, height-7 mm cylinder with `catallax-step::write_step`, then calls `catallax-import::propose`, `apply` and `rebuild` and `catallax_app::build_incremental_with_files`. These entries already existed; this slice adds application-level acceptance evidence, not another importer or a new wire format.

Both sources need no repairs and produce no derivative. Exact carrier, edge-curve and vertex coincidence queries compare the application body with the source primitives. The checks also require sound topology, the file's face origins, matching exact bounding boxes, overlapping volume enclosures narrower than 0.000001 mm³, and a valid checked STEP re-export. Import movement is zero for these fixtures. The source export's reported bound remains a separate quantity and is not inferred from volume agreement.

Display is checked through the application entry at kernel defaults and at SolidDraft's automatic coarse, standard and fine controls, with chord height `0.002 × body diagonal × {4, 1, 0.25}` mm and a 15-degree surface-normal bound. Every face must meet both requested bounds within the existing application work guard. This tests Catallax's controls, not SolidDraft's f32 conversion, rendered viewport or complete application workflow. No outside STEP fixtures are used.

Still refused or unproven by this slice: general cones, pole-bounded spheres and tori, splines, curved gap repair, explicit over-cap acceptance in recipe evaluation, and multiple solids in one recipe import step. The lower-level importer can return separate bodies; the recipe application does not merge them. Seam exchange contracts and field construction from imported faces are unchanged.

SolidDraft uses `propose` for the preview and `apply` for the repair record, saves unchanged source bytes and their SHA-256 plus the recipe decisions, and supplies those bytes through `ImportFiles` to `build_incremental_with_files`. It must use each proposal's exact `question` string for saved decisions, not construct its own gap-index spelling. `rebuild` checks the persisted repair record on reopen.

## Supported first path

Files: AP203/AP214/AP242 `MANIFOLD_SOLID_BREP` with one `CLOSED_SHELL`, in millimetres, inches or metres (exact conversion).

Surfaces: `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`. Curves: `LINE`, `CIRCLE`, `ELLIPSE`, and `SURFACE_CURVE` or `SEAM_CURVE` read through their 3D curve (parameter-space curves are ignored). A `CONICAL_SURFACE` is read only when its semi-angle is exactly 45 degrees in a degree context, because other angles have an irrational slope; others are refused. Several `MANIFOLD_SOLID_BREP` entities import as separate bodies, which acceptance test 2 needs. Everything else (B-spline, offset, swept, voids, shell models, assemblies) is refused with the entity id and type. A face with a single `FACE_BOUND` and no `FACE_OUTER_BOUND` is read with that bound as its outer boundary, as Catallax's own exporter writes it. Decimals read as exact rationals; no float enters a decision.

Face keys are the `ADVANCED_FACE` entity names, as in `docs/RECIPE.md` section 4.18.

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

## Change record (`catallax-import-repair/1`)

The recipe `decisions` saves answers under the proposal's exact `question`, `repair:<sha256-prefix>:g<gap-index>`. The JSON change record is kept beside the recipe so reopening restores it.

```
{
  "format": "catallax-import-repair/1",
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
  }]
}
```

Numbers are decimal strings. `max_movement_mm` is measured on the repaired result, not copied from the proposal. A `user` approval records the approver identity the app supplies and the exact proposed tolerance.

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
- Automatic closure needs a seam, a distance of at most 0.01 mm (inclusive) and a squared distance times 100 not above the squared length of the shortest straight edge of the faces involved. Otherwise the proposal asks, with the reason.
- The repair moves whole vertices of one side onto their partners and refits each planar face through its moved vertices, which must stay exactly coplanar. Both sides are tried. The side that changes fewer surfaces wins, then the smaller movement. So a seam closes by extending faces along their own surfaces where possible. If neither side works, `validation_failed` and nothing is built.
- The repaired solid is rebuilt in a scratch model, checked for closed, sound topology and positive volume before it is built for the caller. Self-intersection is tested exactly (see below) on the source solids and again on every repaired solid.
- Recorded bounds are the bounds of the solid's vertices; they ignore curved bulge. Movement is the largest vertex displacement, rounded up at 1e-12 mm.
- Curved gaps are out of scope for the first path. A gap that involves a curved edge or face is refused with `not_solid` or `validation_failed`, never guessed at. The message names the faces and edges involved, and the gap id once one is assigned. This is a documented limit.
- A decline leaves that solid open and unbuilt: it is reported in `left_open` and not built, the source bytes are never modified, and the application must say plainly that the solid was left open. This is a documented limit.
- A solid with no gap is built as written. Separate solids are never joined, so bodies 0.005 mm apart produce no proposal.

### Self-intersection

`check_self_intersection` and the build entries compare every pair of faces of one solid and refuse with `self_intersecting`, naming the solid and the first pair of faces, or with `self_intersection_not_verified` when a pair cannot be decided. A pass means every pair was decided. Separate solids are never compared with each other.

The rule is the same for every kind of face. Two faces intersect when their carrier surfaces meet in a curve and a positive-length piece of it lies in both closed faces and in the interior of at least one; or when the carriers coincide, the faces look the same way and their interiors overlap. Contact along a shared edge, at a single point (including a sphere tangent to a face) and back-to-back coincident faces facing opposite ways are allowed. A cylinder wall that meets its own cap along their common circle passes. A boss whose ring is a hole boundary of the face it stands on passes. A boss standing on a whole face with no hole for it is refused, because its ring lies in that face's interior.

Planar faces with straight edges are compared exactly in rational arithmetic. Every other pair (cylinder, 45 degree cone, sphere, torus and any face with a circular or elliptical edge) goes to the certified intersection engine, after pairs whose proven boxes are strictly apart are set aside. The engine returns the curves where the carriers meet. Each curve is cut where it meets the edges of either face, one exact point inside every piece is found, and its membership in both faces and on their boundaries is decided with exact or certified answers. No floating-point value decides any step. An engine answer that is not exact is treated as undecided.

Not verified, with the pair and the reason named, and never a pass: a face with no boundary on a plane, cylinder or cone; coincident carriers facing the same way where a region is bounded by curved edges; a boundary edge that lies in both carriers and meets the meeting curve without matching it; any pair the engine refuses or cannot make exact; and a solid whose faces cannot be set up for the engine. An intersection found for any pair wins over an undecided pair. A solid with an undecided pair is refused, so it cannot be imported until the pair is decidable.

Not fixtured: a repair that itself creates an intersection. The repaired path shares the checker, and no fixture forces such a repair.
