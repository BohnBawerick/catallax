# Capabilities and limits

This page describes the Catallax v0.5.0 and SolidDraft v0.8.0 cut. Use the [roadmap](../ROADMAP.md) for work that is still planned. A supported operation can still refuse outside its checked domain. See the [recipe](RECIPE.md) for those domains and the [verification protocol](VERIFY-INTERFACE.md) for what an answer proves.

## What works today

- The B-rep builds stock solids, sketches, named parametric features and Booleans in supported configurations. It stores analytic and spline definitions and keeps operation-derived face and edge names through supported edits. Shell, draft, variable fillets, freeform sweeps and spline cuts still have gaps.
- STEP import distinguishes admitted solids, partial openings and retained surfaces. Unproved surface previews do not count as editable solids. Checked STEP export writes analytic geometry or labels a fitted result approximate. Assemblies can carry shared parts, placements, names, colours and supported semantic PMI. IGES, closed-mesh import and scan fitting cover limited domains.
- Geometry queries cover supported mass, distance, clash, clearance, section, silhouette and hidden-line cases. General cone and torus intersections and all three phase 2 target chains still refuse in the release scoreboard.
- Both evaluators use the same recipe definitions. The field builds supported lattices and gyroids, with separate print and mass stand-in requests. [The comparison route](COMPARE.md) reports sampled agreement, disagreement or unresolved checks. It does not certify a whole body.
- SolidDraft uses Catallax for modelling, assemblies and drawings. A kernel library result does not by itself establish app acceptance. The [CLI](BUILDING.md) can replay requests and compare strict output bytes across worker counts.

## Release evidence

| Measure | Catallax | OpenCASCADE | SMLib |
| --- | --- | --- | --- |
| Edit chains, P0 to P5 | 11 of 13 | 11 of 13 | 7 of 13 |
| Edit chains, phase 2 targets | 0 of 3 | 3 of 3 | 3 of 3 |
| Edit chains, census parts | 9 of 10 | 10 of 10 | 8 of 10 |
| STEP census, complete-source valid solid exports | 767 of 3,000 | 2,701 healed-valid reference | No run |

Full solid admission is a separate count: 889 of 3,000. The older OpenCASCADE reference used healing and a different time limit. It does not establish a same-protocol import gap. The [coverage reference](COVERAGE.md) defines the 61 lines; the release scored 31 points, with one full line and 60 partial lines.

## Imported-solid direct edits

`catallax_import::direct::DirectSolid` holds the admitted source body. `apply` builds a private candidate and publishes it only after sewing and numerical admission. A refusal leaves the source and edit history unchanged. Use `body()` for export, not `source_faces()`.

| Request | Supported construction |
| --- | --- |
| `DeleteFace`, `DeleteFillet` | Remove a four-sided face or cylindrical blend between planar side supports and parallel end planes. |
| `DeleteHole` | Remove connected inward coaxial cylindrical walls and blind or stepped caps, then remove complete mouth loops. |
| `MoveFace`, `OffsetFace`, `ReplaceFace` | Rebuild a planar or cylindrical support and its shared boundaries. Cylinder caps must be perpendicular to the axis. |

Edits retain face correspondence and record sewing and endpoint corrections separately from the requested movement. A rational construction does not make a numerical imported body exact. General spline replacement and higher-valence corner reconstruction are outside these paths. These are library APIs, not SolidDraft commands.

## Display meshes

`catallax_mesh::mesh_body` builds a mesh tagged with input model face IDs. The default chord tolerance is 0.01 mm. Faces share one polyline per edge. The report records chord error and available angular, edge-length and aspect-ratio checks. A resource or precision limit refuses the mesh rather than returning an unchecked result at a looser tolerance. `mesh_body_reusing` can reuse a face when the supporting definition, trim, boundary samples and controls match. See the [recipe's application entry](RECIPE.md#63-application-entry) for picking and revision details.
