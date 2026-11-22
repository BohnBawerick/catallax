# Coverage against OpenCASCADE and SMLib

The combined scope is 61 capability lines. Full is 1, partial is 0.5, and not started is 0. A general capability must work for full credit; a limited domain or library call with real refusals is partial. Catallax-only capabilities, including the independent field and proof button, add no points. Speed does not fail a change before phase 5.

## Sources and limits

Adopted by [Checkpoint 2 decisions, Q2 and section 4.10](checkpoints/checkpoint-2-decisions.md#410-coverage-the-12-changes-for-100). The list and historical scores come from the Checkpoint 2 review report, Appendix A, "The 61 lines", also reproduced in the grill pack Part 3. Full requirements below come from that report's "Smallest set of additions for 100% of both kernels", the grill pack Part 3, and their existing roadmap done-when lines. Those private working documents are not bundled; their adopted twelve changes are preserved in the [roadmap coverage follow-on](../ROADMAP.md#coverage-follow-on-complete-the-61-lines). [The checkpoint report](checkpoints/checkpoint-2.md) holds the measured release results.

The line list was reconstructed on 19 Nov 2026 at `ade010bf` to fit the published area totals. It is not a recovered original 14 Nov list. Reference-kernel membership and the per-line 14 Nov scores were inferred in that review. Historical scores below are not new verification and do not certify the new full bars.

A full bar that the sources do not state is **bar to be set**, even for a line previously scored full. Known requirements are retained beside that marker, not promoted into an invented complete test. Such a line cannot be newly claimed full until its bar is set. Generated sets and census acceptance must be frozen before implementation; this document does not invent sample counts or pass rates.

Kernels: B = both OpenCASCADE and SMLib; O = OpenCASCADE only; S = SMLib only. "19 Nov" is the historical review score, not current acceptance. Each roadmap done-when must name its coverage IDs and level; the twelve coverage extensions target full, 1 point.

## Curves and surfaces (12)

| ID | Capability | Kernels | 14 Nov | 19 Nov | Full bar |
|---|---|---|---:|---:|---|
| C1 | Exact lines, circles, ellipses and other conics | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| C2 | NURBS curves and surfaces with exact rational evaluation | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| C3 | Analytic surfaces: plane, cylinder, cone, sphere, torus, elliptic cylinder | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| C4 | 2D sketch geometry: lines, arcs, tangent arcs, trims, 2D fillets | O | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| C5 | Edit a spline: insert knots, raise degree, split, join | B | 0.5 | 1 | Knot insertion, degree raise, split and join are public for spline editing. Source: existing 2.2 done-when. |
| C6 | Offset curves and offset surfaces | B | 0.5 | 0.5 | Offset any curve or surface, spline faces included, within a stated tolerance. Change 2. |
| C7 | Fit a spline curve through or near points | B | 0.5 | 0.5 | Least-squares fitting of curves to points within a stated tolerance; each builds on a generated set without refusing possible geometry. Change 1. |
| C8 | Surface through a grid of points, or fitted to a point cloud | B | 0 | 0.5 | Least-squares fitting of surfaces to points within a stated tolerance, alongside grid interpolation; each builds on a generated set without refusing possible geometry. Change 1. |
| C9 | Surfaces from curve networks: Coons, Gordon, guided skin, N-sided fill | B | 0 | 0.5 | G1 N-sided and Gordon fills, and curve networks whose corners miss by up to the stated tolerance; each builds on a generated set without refusing possible geometry. Change 1. |
| C10 | Skinned and swept surfaces through sections or along a path | B | 0.5 | 0.5 | bar to be set. Stated extension: lofts with guide curves build on a generated set without refusing possible geometry, change 1. A complete bar for all skinned and swept surfaces is not stated. |
| C11 | Trim, extend and project curves and surfaces | B | 0.5 | 0.5 | bar to be set. Stated extension: project curves onto surfaces, change 2. The sources do not state the complete general trim-and-extend bar for curves and surfaces. |
| C12 | Bend and twist a body | S | 0 | 0 | Bend or twist any body, general input rather than a first build. Change 3. |

## Intersection (6)

| ID | Capability | Kernels | 14 Nov | 19 Nov | Full bar |
|---|---|---|---:|---:|---|
| I1 | Analytic surface pairs: plane, cylinder, sphere, cone, torus | B | 0.5 | 0.5 | Every analytic pair, including skew cone/torus and cylinder/torus; no refusal on the generated set. Change 4. |
| I2 | Spline surface against analytic surface | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| I3 | Spline surface against spline surface | B | 0.5 | 0.5 | Spline against spline with no refusal on the generated set. Change 4. |
| I4 | Tangent and near-tangent contact | B | 0.5 | 0.5 | Tangent and near-tangent contact answered under the bounded-answer rule. Change 4. |
| I5 | Curve against surface and curve against curve | B | 0.5 | 0.5 | bar to be set. Curve against surface must have no refusal on the generated set, change 4; the complete curve-against-curve bar is not stated. |
| I6 | Plane and surface sections of a solid | B | 0.5 | 0.5 | bar to be set. Existing 3.1 requires exact section curves, silhouettes and hidden-line views used by SolidDraft instead of mesh outlines, but no separate full I6 bar is stated. |

## Solids and topology (7)

| ID | Capability | Kernels | 14 Nov | 19 Nov | Full bar |
|---|---|---|---:|---:|---|
| S1 | Primitive solids: block, cylinder, sphere, cone, torus | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| S2 | Unite, subtract and intersect analytic solids | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| S3 | Booleans on solids with spline faces | B | 0.5 | 0.5 | Booleans on spline-faced solids pass a generated set and the census parts, not just named requests. Change 5. No clock or 20M gate before phase 5. |
| S4 | Several bodies in one part: combine, move, copy | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| S5 | Split a solid into pieces with a plane or tool body | O | 0.5 | 0.5 | Split a solid by a plane, sheet or body, keeping every piece. Change 6. |
| S6 | Sew faces into a closed solid and check it | B | 0.5 | 0.5 | Sew any face set that closes into a checked solid. Change 6. |
| S7 | Surface, wire and non-manifold bodies | B | 0 | 0.5 | Surface, wire and non-manifold bodies. General surface-body scope, not sheets and wires only. Change 6. |

## Features that add or remove material (11)

| ID | Capability | Kernels | 14 Nov | 19 Nov | Full bar |
|---|---|---|---:|---:|---|
| F1 | Extrude a profile: blind, two-sided, through all | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| F2 | Revolve a profile, full or partial | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| F3 | Linear and circular patterns, mirror | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| F4 | Extrude up to a face or to the next face | O | 0.5 | 0.5 | Up-to-face builds on line, arc, ellipse and spline profiles, up to any curved face. Change 7. |
| F5 | Extrude with taper | B | 0 | 0.5 | Taper builds on line, arc, ellipse and spline profiles. Change 7. |
| F6 | Sweep a profile along a path | B | 0.5 | 0.5 | Sweep any line, arc, ellipse or spline profile along any path. Change 7. |
| F7 | Sweep with a guide, with twist or along a helix | B | 0 | 0.5 | Guide, twist and helix each build on line, arc, ellipse and spline profiles along any path. Change 7. P8 remains parked, not silently counted done. |
| F8 | Loft through sections | B | 0.5 | 0.5 | Lofts with guide curves build on a generated set without refusing possible geometry. Change 1. |
| F9 | Pipe: a section along a 3D curve | B | 0.5 | 0.5 | A section along any 3D path, including guide, twist and helix, with line, arc, ellipse and spline profiles. Change 7. |
| F10 | Imprint curves on faces and split faces | B | 0.5 | 0.5 | Imprint a curve on a face and split faces as a user operation. Change 6. |
| F11 | Extrude, revolve and sweep spline and ellipse profiles | B | 0.5 | 0.5 | bar to be set. Change 7 requires the general taper/sweep/up-to-face family on line, arc, ellipse and spline profiles; a complete full bar including revolve is not separately stated. |

## Edge and wall features (9)

| ID | Capability | Kernels | 14 Nov | 19 Nov | Full bar |
|---|---|---|---:|---:|---|
| E1 | Fillet with constant radius on any edge | B | 0.5 | 0.5 | Constant-radius fillets on any edge; the census parts take the feature without refusal of possible geometry. Change 8. |
| E2 | Fillet with variable radius | B | 0.5 | 0.5 | Variable-radius fillets with any law on any edge; the census parts take the feature without refusal of possible geometry. Change 8. |
| E3 | Fillet between spline faces, across patches | B | 0.5 | 0.5 | Fillets on spline-faced solids, across patches, pass a generated set and the census parts. Change 5. Existing 2.1 still requires the radius edit to keep names and the curved-housing chain with later edits and unhealed OCCT STEP. |
| E4 | Chamfer with equal or unequal setbacks | B | 0.5 | 0.5 | Chamfer roll-over and general corners, with equal or unequal setbacks; the census parts take the feature without refusal of possible geometry. Change 8. |
| E5 | Shell: hollow a solid with open faces | B | 0.5 | 0.5 | Shell on any face, spline and oblique cylinder included; the census parts take the feature without refusal of possible geometry. Change 8. |
| E6 | Offset or thicken a solid or a sheet | B | 0.5 | 0.5 | Offset any whole solid or surface, spline faces included, within a stated tolerance; 2.3 also requires thickening a surface into a solid. Changes 2 and 6. |
| E7 | Draft faces | O | 0.5 | 0.5 | Draft on any face, spline and oblique cylinder included; the census parts take the feature without refusal of possible geometry. Change 8. |
| E8 | Curvature-continuous (G2) blends | S | 0 | 0 | Curvature-continuous blends on general input, not a first build. Change 3. |
| E9 | Remove a face, fillet or hole; move, offset or replace a face | O | 0 | 0.5 | Each direct edit works on any face type, through SolidDraft. Change 9. |

## File exchange (8)

| ID | Capability | Kernels | 14 Nov | 19 Nov | Full bar |
|---|---|---|---:|---:|---|
| X1 | STEP export of solids | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| X2 | STEP import of real solids | B | 0.5 | 0.5 | bar to be set. Change 10 requires STEP import parity against OpenCASCADE's 90.03% healed-valid census reference; the exact census target and comparable protocol must be set, not invented. |
| X3 | IGES import and export | B | 0.5 | 0.5 | bar to be set. Existing 3.2 adds IGES sphere and torus export; the complete IGES import/export bar is not stated. |
| X4 | Repair imported geometry: gaps, tolerances, sewing | B | 0.5 | 0.5 | bar to be set. Change 10 requires STEP import parity against OpenCASCADE's 90.03% healed-valid reference. The complete repair acceptance and census target are not stated. |
| X5 | Mesh export: STL, OBJ, 3MF | O | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| X6 | Mesh import into a solid | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| X7 | Surface, wire and multi-body files | B | 0.5 | 0.5 | Surface, wire and multi-body exchange; sheets and wires export as sewn shells, not loose geometric sets. Change 6. |
| X8 | Assemblies, colours, names and annotations in files | O | 0 | 0 | Assemblies, colours, names and annotations round-trip in STEP. Change 11. The annotations include the P1 semantic-PMI C1 to C5 acceptance; drawings alone do not close the rest of X8. |

[Item 3.2](milestones/3.2-assemblies-and-exchange.md) records the kernel exchange checks and remaining scope for X3 and X8.

## Meshes and queries (8)

| ID | Capability | Kernels | 14 Nov | 19 Nov | Full bar |
|---|---|---|---:|---:|---|
| Q1 | Display meshes with chord and angle tolerances | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| Q2 | Mass properties: volume, area, centroid, inertia | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| Q3 | Point inside or outside a solid | B | 1 | 1 | bar to be set. The sources name the capability and historical score but do not state a complete acceptance bar. |
| Q4 | Closest point on a curve or surface | B | 0.5 | 0.5 | Closest point on any curve or surface. Change 12. |
| Q5 | Bounds, curvature and local surface analysis | B | 0.5 | 0.5 | Bounds, curvature and local surface analysis on general input. Change 12. |
| Q6 | Clash, clearance and minimum distance between solids | O | 0.5 | 0.5 | No pair answers "unresolved" outside the stated tolerance, including distant and nested pairs. Change 12 widens 3.3. |
| Q7 | Hidden-line views | O | 0 | 0 | SolidDraft uses exact kernel section curves, silhouettes and hidden-line views in place of mesh outlines. Existing 3.1. |
| Q8 | Mesh Booleans and decimation | S | 0 | 0 | Mesh Booleans and decimation on general input, not a first build. Change 3. |

## Historical totals, not a new recount

| Measure | 14 Nov, reconstructed | 19 Nov review |
|---|---:|---:|
| Points of 61 | 33.5 | 37.0 |
| Combined scope | 55% | 61% |
| Of OpenCASCADE's 58 lines | 58% | 64% |
| Of SMLib's 52 lines | 57% | 62% |
| Complete / partial / not started | 17 / 33 / 11 | 18 / 38 / 5 |

The review estimated about 67% if the old milestone demonstrations all passed, and 85% at their full planned scope. Those were planning estimates, not evidence of 100%. The twelve approved changes target full coverage; unknown bars remain explicit. Rule R5 permits one scoreboard, census and coverage recount per checkpoint, no dry runs.
