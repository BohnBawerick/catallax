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
| Declined | Original body used unchanged; no derivative is made; the decision is saved as `decline` so the question is not re-asked for the same bytes |

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
