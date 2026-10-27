# Piece 8: STEP import and bounded gap repair (contract)

Status: contract only. No importer exists on main. `docs/DESIGN.md` records "The importer is not built yet", `catallax-step` only writes, and `Operation::ImportStep` in `catallax-recipe` is a schema entry with no consumer. This file fixes the narrow first supported path so the build can start from written text. It does not claim any implemented behavior.

## Why this is not a slice yet

A real click path needs, in order: a Part 21 reader, analytic surface and curve recognition, exact unit conversion, topology construction with periodic faces and sharp edges, gap measurement, a fixed imported body in the recipe consumer, field construction from the faces, the persisted change record, and the SolidDraft dialog. The topology, recipe consumer and app adapter files belong to other lanes (L01, L03, L04, L07). A reader plus a dialog without the body path would be a stub, so none is shipped.

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

Stored in the recipe `decisions` under the question id `repair:<sha256-prefix>:<gap-index>`, with the JSON record kept beside the recipe so reopening restores it.

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

`check_self_intersection` and the build entries compare every pair of planar, straight-edged faces of one solid in rational arithmetic and refuse with `self_intersecting`, naming the solid and the first pair of faces. Two faces intersect when their planes meet in a line and a positive-length piece of it lies in both closed faces and in the interior of at least one; or when they are coplanar, face the same way and their interiors overlap. Contact along a shared edge, at a single vertex or edge point, and back-to-back coplanar faces facing opposite ways are allowed. `apply` checks the source solids first, including ones that stay open, and the scratch build of each repaired solid checks it again, returning the typed refusal itself.

Limits: pairs involving a curved surface or a circular edge are not compared, so a cylinder cap or wall is not tested against other faces. Separate solids are never compared with each other. The repaired-solid path shares the same checker, but no fixture forces a repair that creates an intersection.
