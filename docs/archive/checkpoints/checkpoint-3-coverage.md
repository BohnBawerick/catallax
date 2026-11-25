Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Coverage](../../COVERAGE.md).

# Checkpoint 3 coverage recount

This recount scores the [61 existing lines](../../COVERAGE.md) at [Catallax v0.5.0](https://github.com/BohnBawerick/catallax/releases/tag/v0.5.0), with application evidence from [SolidDraft v0.8.0](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/releases/tag/v0.8.0), which pins the same kernel. Each line earns 1 for its documented full bar, 0.5 for working evidence within a limited domain, or 0 for no working evidence. The release receipts below support the scores. [The checkpoint scoreboard](checkpoint-3.md#edit-chain-scoreboard) and [STEP census](checkpoint-4.md#step-census) give part and file results separately.

Sixteen lines that scored 1 in the 19 Nov review still have no complete written bar. They earn 0.5 in this recount. Q1 has a written bar and also earns 0.5 because mesher acceptance remains open. Bend and twist, G2 fillets, assembly exchange, hidden-line views and mesh Booleans and decimation gain half a point each. The new total uses the adopted acceptance bars, so a drop from the historical score alone does not establish a capability regression.

| Scope | Points | Share |
|---|---:|---:|
| Combined | 31 of 61 | 50.8% |
| OpenCASCADE | 29.5 of 58 | 50.9% |
| SMLib | 26.5 of 52 | 51.0% |

One line meets its full bar, 60 have partial evidence and every line has working evidence. The independent field and the phase-5 work targets add no points.

## Curves and surfaces

| ID | Score | Evidence at the cut |
|---|---:|---|
| C1 | 0.5 | Exact conics are represented and evaluated. The full bar still needs a definition. [Kernel operation domains](../CAPABILITIES-before-checkpoint-3.md#what-works-today) |
| C2 | 0.5 | Rational spline definitions have exact evaluation. The full bar still needs a definition. [Kernel operation domains](../CAPABILITIES-before-checkpoint-3.md#what-works-today) |
| C3 | 0.5 | Analytic carriers include elliptic cylinders. The full bar still needs a definition. [Kernel operation domains](../CAPABILITIES-before-checkpoint-3.md#what-works-today) |
| C4 | 0.5 | Solved sketches carry lines, arcs and construction constraints. The full bar still needs a definition. [Sketch contract](../../RECIPE.md#21-sketch) |
| C5 | 1 | Knot insertion, degree raise, split and join are public. [NURBS toolkit](https://github.com/BohnBawerick/catallax/pull/676) |
| C6 | 0.5 | Supported offsets include planar and cylindrical faces and signed spline-sheet construction. General curve and surface offsets remain open. [Offset contract](../../RECIPE.md#48-offset) [Surface operations](../../SURFACE-BODIES.md#public-operations) |
| C7 | 0.5 | Spline interpolation works. The generated-set requirement for least-squares fitting remains open. [NURBS toolkit](https://github.com/BohnBawerick/catallax/pull/676) |
| C8 | 0.5 | Point-grid interpolation works. The sphere scan fit does not cover general spline fitting from a point cloud. [NURBS toolkit](https://github.com/BohnBawerick/catallax/pull/676) [Scan result](../milestones/4.4-imports-and-scans.md#measured-scan) |
| C9 | 0.5 | Coons, Gordon, guided skin and N-sided fills have working examples. The N-sided fill is C0 and the generated-set and corner-gap bars remain open. [NURBS toolkit](https://github.com/BohnBawerick/catallax/pull/676) |
| C10 | 0.5 | Supported lofts and path sweeps build. Guided loft coverage on the generated set remains open. The full bar still needs a definition. [Loft](../../RECIPE.md#412-loft) [Sweep](../../RECIPE.md#411-sweep) |
| C11 | 0.5 | Rectangular surface trim and tangent extension work. General curve projection and a complete trim-and-extend bar remain open. [Surface operations](../../SURFACE-BODIES.md#public-operations) |
| C12 | 0.5 | Bend and twist each build on supported prismatic inputs. General bodies remain open. [Bend and twist contract](../../RECIPE.md#420-bend-and-twist) |

## Intersection

| ID | Score | Evidence at the cut |
|---|---:|---|
| I1 | 0.5 | Plane, cylinder and sphere pairs and some cone and torus cases work. The scoreboard's general cone and torus rows still refuse. [Kernel operation domains](../CAPABILITIES-before-checkpoint-3.md#what-works-today) [Scoreboard](checkpoint-3.md#edit-chain-scoreboard) |
| I2 | 0.5 | Spline faces meet supported analytic cutters. The full bar still needs a definition. [Spline cuts](../phase1/spline-boolean-top-cuts/REPORT.md) |
| I3 | 0.5 | Certified spline-pair branches work on supported pairs. Tangent or unseparated branches still refuse. [Spline intersection research](../research/spline-intersection-2026-method.md) |
| I4 | 0.5 | Typed tangent results and tolerance-aware crossing checks exist. General bounded answers near contact remain open. [P0 controls](checkpoint-3.md#edit-chain-scoreboard) [Tolerance-aware pair checks](https://github.com/BohnBawerick/catallax/pull/741) |
| I5 | 0.5 | The pair search handles supported straight and circular edges and crossings. The generated curve-surface set and complete curve-curve bar remain open. [Pair search](https://github.com/BohnBawerick/catallax/pull/741) |
| I6 | 0.5 | Face-by-face planar sections and drawing curves ship, with exact and fitted answers. General surface sections have no complete bar. [Section and drawing curves](../interface-changes/2026-11-21-sections-and-drawing-curves.md) |

## Solids and topology

| ID | Score | Evidence at the cut |
|---|---:|---|
| S1 | 0.5 | Block, cylinder, sphere, cone and torus stock solids exist. The full bar still needs a definition. [Stock solids](../../RECIPE.md#40-stock-solids) |
| S2 | 0.5 | Analytic unite, subtract and intersect work on supported pairs. The full bar still needs a definition. [Boolean contract](../../RECIPE.md#44-boolean) |
| S3 | 0.5 | The recorded loft-top cuts reach checked approximate STEP, and one through-hole family continues after a fitted fillet. The generated set and census-part acceptance remain open. [Loft cuts](../phase1/spline-boolean-top-cuts/REPORT.md) [One model](../milestones/3.0B-one-model.md) |
| S4 | 0.5 | Multiple placed bodies and a rooted assembly work. The full bar still needs a definition. [Recipe placements](../../RECIPE.md#13-numbers-units-and-placements) [Assembly exchange](../milestones/3.2-assemblies-and-exchange.md) |
| S5 | 0.5 | Plane sections and supported half-space cuts work. General solid split by a sheet or tool with every piece returned remains open. [Sections](../interface-changes/2026-11-21-sections-and-drawing-curves.md) [Boolean contract](../../RECIPE.md#44-boolean) |
| S6 | 0.5 | The six-face sewn box and retained-import sewing have working evidence. General open sets that can close still lack acceptance. [Sewing proof](../../SURFACE-BODIES.md#reproduction-and-observed-proof) |
| S7 | 0.5 | Sheets and free wires have kernel and exchange representations. General non-manifold bodies remain open. [Body roots](../../SURFACE-BODIES.md#body-roots-and-step) |

## Features that add or remove material

| ID | Score | Evidence at the cut |
|---|---:|---|
| F1 | 0.5 | Supported profile extrusions include blind, two-sided and through-all variants. The full bar still needs a definition. [Extrusion contract](../../RECIPE.md#41-extrude) |
| F2 | 0.5 | Full and partial revolves have a recipe operation. The full bar still needs a definition. [Revolve contract](../../RECIPE.md#42-revolve) |
| F3 | 0.5 | Linear and circular patterns and mirror have recipe operations. The full bar still needs a definition. [Pattern](../../RECIPE.md#413-pattern) [Mirror](../../RECIPE.md#414-mirror) |
| F4 | 0.5 | Recipe 6 resolves supported trimmed planar and curved start and end faces. General curved limits across all profile types remain open. [Face-bounded extrusion](../interface-changes/2026-11-24-face-bounded-extrusion.md) |
| F5 | 0.5 | Tapered extrusion works within the supported domain. The complete line, arc, ellipse and spline matrix remains open. [Tapered extrusion](../../RECIPE.md#413-tapered-extrusion) |
| F6 | 0.5 | Sweeps along supported lines, arcs and splines work. The all-profiles, all-paths bar remains open. [Sweep contract](../../RECIPE.md#411-sweep) |
| F7 | 0.5 | Guided, twisted and helical sweeps have working cases. The general profile and path matrix remains open, and P8 still refuses before the spring. [Numerical sweeps](../../RECIPE.md#4112-numerical-guides-twist-and-helices) [Scoreboard](checkpoint-3.md#edit-chain-scoreboard) |
| F8 | 0.5 | Lofts through supported sections build. The generated guided-loft set remains open. [Loft contract](../../RECIPE.md#412-loft) |
| F9 | 0.5 | Supported 3D sweep paths include helical builds. Arbitrary paths and profiles lack acceptance. [Numerical sweeps](../../RECIPE.md#4112-numerical-guides-twist-and-helices) |
| F10 | 0.5 | The Boolean builder imprints and splits faces internally. A general user operation remains open. [Kernel operation domains](../CAPABILITIES-before-checkpoint-3.md#what-works-today) |
| F11 | 0.5 | Selected spline and ellipse profiles work in the feature builders. General sweep coverage and a complete revolve bar remain open. [Recipe profiles and paths](../../RECIPE.md#2-profiles-and-paths) |

## Edge and wall features

| ID | Score | Evidence at the cut |
|---|---:|---|
| E1 | 0.5 | Supported constant-radius fillets include strips, rims and selected complex contacts. The census-wide any-edge bar remains open. [Fillet contract](../../RECIPE.md#45-fillet) |
| E2 | 0.5 | Variable-radius fillets work with supported laws and carriers. General laws and census edges remain open. [Variable-radius fillets](../../RECIPE.md#452-variable-radius-fillets) |
| E3 | 0.5 | Multi-patch spline fillets build in selected cases. Imported coverage, edited names and the whole curved-housing chain remain open. [Checkpoint 2 evidence](checkpoint-2.md#what-phase-2-delivered-against-roadmapmd) [Scoreboard](checkpoint-3.md#edit-chain-scoreboard) |
| E4 | 0.5 | Supported chamfers have equal and unequal setbacks. General roll-over, corners and census acceptance remain open. [Chamfer contract](../../RECIPE.md#46-chamfer) |
| E5 | 0.5 | The shelled tray builds. General spline and oblique-cylinder shell acceptance remains open. [Shell contract](../../RECIPE.md#47-shell) [Scoreboard](checkpoint-3.md#edit-chain-scoreboard) |
| E6 | 0.5 | Selected offsets and signed thickening produce solids. Offset on any spline face and general thickening remain open. [Thickening](../../RECIPE.md#49-thicken) [Surface operations](../../SURFACE-BODIES.md#public-operations) |
| E7 | 0.5 | Planar and selected cylindrical draft cases work. General face and census coverage remain open. [Draft contract](../../RECIPE.md#410-draft) |
| E8 | 0.5 | G2 fillets have a first build in a limited support domain. The general-input bar remains open. [G2 fillets](../../RECIPE.md#453-curvature-continuous-fillets) |
| E9 | 0.5 | The library has the six direct edits for supported neighbourhoods. Any-face support through SolidDraft remains open. [Direct edit limits](../CAPABILITIES-before-checkpoint-3.md#imported-solid-direct-edits) |

## File exchange

| ID | Score | Evidence at the cut |
|---|---:|---|
| X1 | 0.5 | Analytic and checked fitted solids write STEP. The full bar still needs a definition. [Kernel operation domains](../CAPABILITIES-before-checkpoint-3.md#what-works-today) [Scoreboard](checkpoint-3.md#edit-chain-scoreboard) |
| X2 | 0.5 | The frozen census measures full and partial solid admissions and valid exported files. The parity target and comparable protocol still need a definition. [Census](checkpoint-4.md#step-census) |
| X3 | 0.5 | IGES reads supported geometry and writes whole spheres and ring tori. General import/export acceptance remains open. [IGES import](../../RECIPE.md#4182-iges-import) [Exchange checks](../milestones/3.2-assemblies-and-exchange.md#coverage) |
| X4 | 0.5 | The gap rule and per-body repair admit some checked and tolerant bodies. General repair acceptance and the parity target remain open. [Import gap rule](../milestones/3.0A-import-gap-rule.md) [Census](checkpoint-4.md#step-census) |
| X5 | 0.5 | STL, OBJ and 3MF export paths exist. The full bar still needs a definition. [SolidDraft release](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/releases/tag/v0.8.0) [Printable lattice](../milestones/4.3-printing.md#p6-at-001-mm) |
| X6 | 0.5 | Closed mesh inputs can become solids. The full bar still needs a definition. [Mesh import](../../RECIPE.md#419-mesh-import) |
| X7 | 0.5 | STEP carries supported sheets, wires and multiple bodies. General sewn-shell exchange remains open. [Body roots and STEP](../../SURFACE-BODIES.md#body-roots-and-step) [Assembly exchange](../milestones/3.2-assemblies-and-exchange.md) |
| X8 | 0.5 | A rooted STEP assembly round-trips with shared definitions, per-part colours, names and supported PMI. Wider encodings, occurrence and face styles, and complete P1 PMI reader and browser checks remain open. [Exchange checks](../milestones/3.2-assemblies-and-exchange.md#coverage) [PMI exits](checkpoint-3.md#phase-3-delivery) |

## Meshes and queries

| ID | Score | Evidence at the cut |
|---|---:|---|
| Q1 | 0.5 | Shared charts produce bounded meshes in supported cases. The full frozen corpus, shared-edge checks, dense bounds, point prediction and cross-machine proof remain open. [Mesher audit](checkpoint-3.md#one-display-mesher) |
| Q2 | 0.5 | Volume, area, centroid and inertia have checked kernel paths. The full bar still needs a definition. [Kernel operation domains](../CAPABILITIES-before-checkpoint-3.md#what-works-today) |
| Q3 | 0.5 | Solid membership works on supported families. The full bar still needs a definition. [Scoreboard](checkpoint-3.md#edit-chain-scoreboard) |
| Q4 | 0.5 | Closest-point queries cover supported definitions and trimmed domains. Arbitrary curves and surfaces remain open. [Nearest named face](../../FIELD.md#nearest-named-face) |
| Q5 | 0.5 | Bounds, normal and curvature queries work on supported charts. General local surface analysis remains open. [Kernel operation domains](../CAPABILITIES-before-checkpoint-3.md#what-works-today) [Display bounds](../CAPABILITIES-before-checkpoint-3.md#display-meshes) |
| Q6 | 0.5 | The test assembly gets 38 clear, two clash and 38 unresolved answers. The full bar requires resolving pairs outside tolerance. [Assembly results](../milestones/3.3-clash-clearance.md#results-3) |
| Q7 | 0.5 | SolidDraft drawings use kernel curves. Some curves are fitted and labelled approximate, and unsupported views can refuse. The general exact-view bar remains open. [SolidDraft drawing acceptance](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/357) |
| Q8 | 0.5 | Mesh Booleans and decimation have first builds and explicit stop conditions. General-input acceptance remains open. [Mesh Booleans](../../RECIPE.md#421-mesh-booleans) [Decimation](../../RECIPE.md#422-mesh-decimation) |
