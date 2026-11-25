Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Role of the Catallax field (SDF) kernel under a tolerance-first B-rep

Worktree: Catallax repository at bb2f8375 (abbreviated CX below).
SolidDraft: SolidDraft repository at 0098b8eb (0.7.1 alpha on Catallax v0.4.0) (SD below).
Grill pack: Checkpoint 2 grill pack, private working source (GP below).
Labels: [read] = I read the file or page text; [snippet] = search-result summary only; [unread/memory] = from my own knowledge, not checked this session.

## 1. What the field evaluator does today

### Input and model
- Compiles the recipe, never the B-rep. "It calls no B-rep builder" (CX docs/FIELD.md:5) [read]. Fillets: "The field resolves the source supports from recipe origins ... It constructs no B-rep" (FIELD.md:348) [read]. This confirms the earlier note: it works on recipe definitions, not on B-rep edges.
- Supported (FIELD.md:5, CAPABILITIES.md:36) [read]: stock solids, line/circle/arc sketch profiles, signed extrusions, through-all cuts, full and partial revolutions, named holes (incl. cones), planar face/next extrusions, linear and circular patterns, mirrors, regularised union/subtract/intersect, planar and cylindrical draft (FIELD.md:180-290), elementary fillets (plane-plane cylinder strips, cylinder-plane, plane-perpendicular-cylinder torus strips, three-plane corners; FIELD.md:346-352), chamfers (CAPABILITIES.md:41), offset/shell/thicken analytic reductions (FIELD.md:529), reduced sweeps and lofts with conservative residuals (FIELD.md:310-336, M18-NATIVE-FIELDS.md), closed-mesh import with exact ray parity (FIELD.md:11), STEP import limited to analytic faces, sphere/torus, flat bilinear patches and saddle quadrics; "wider curved trims and source families retain typed coverage refusals" (FIELD.md:~140) [read].
- Field-only features: periodic beam lattices (cubic, constant or linear radius grade), network gyroids (constant or linear level along one axis), explicit positive-level sheet gyroids for queries only (FIELD.md:9-21, 43-46) [read].
- Refused: curved-curved fillets (saddle, oblique, coaxial, fillet-over-fillet) with `field.fillet-contact`; P5 radius-2 fillet with `fillet.unsupported-overflow` (M18-P5-VERDICT.md:13-25) [read]; distance grades, reflected frames, non-axis grades (FIELD.md:13, 17); disconnected gyroid clips `field.print-disconnected` (FIELD.md:84); printing of anything outside primitives/analytic cells/convex planar extrusions/beams `field.print-coverage` (FIELD.md:25); "General fillets, FEA and application integration are not implemented" (FIELD.md:7) [read].

### Outputs
- Point queries returning `Value::TrueDistance` (certified Euclidean distance enclosure) or `Value::Residual` (not a distance, for sweeps and Boolean min/max) (FIELD.md:152-153) [read]. So it is a proper SDF only on primitives and analytic drafts; after Booleans and sweeps it is a conservative implicit, as the PRD admits (GP:3895).
- Exact inside/outside/boundary classification (FIELD.md:402-412).
- Box classifier Inside/Outside/Unknown (FIELD.md:41-58).
- Certified watertight STL/OBJ at a stated two-sided distance bound, binary32-exact coordinates (FIELD.md:23-39, 60-90) [read].
- Proven mass, centroid, inertia enclosures of the field definition (FIELD.md:92-121); mass stand-in STEP of proxy boxes (FIELD.md:106-115).
- Unproved preview meshes, labelled `Proof::Unproved` (FIELD.md:90; CX crates/catallax-app/src/preview.rs:1,10) [read].
- Trimmed minimum distances, nearest named face, clearance (FIELD.md:542-781; docs/M18-FINITE-CLEARANCE.md) [read headings only].
- Scan fitting `reconstruction::scan::fit_scan` lives in this crate but is least-squares point-cloud fitting, not SDF work; the result is a recipe built by both kernels (docs/M23-SCAN-FIT.md:1-20) [read].

### How it is used
- Measuring stick, in tests: catallax-features has catallax-field as a dev-dependency only (crates/catallax-features/Cargo.toml:36-38) [read]. 32 feature test files import it (my grep). Rough tally (grep for Location::Inside/Outside and field refusal codes, heuristic): agreement checks dominate in draft, thicken, holes, cone/cylinder distance, revolutions, open spherical cap (e.g. draft_preserved_boss 44 hits, cone_distance 45); the curved-fillet files (saddle_fillets, oblique_fillets, coaxial_fillets, fillet_over_fillet, offset_fit) only assert that the field refuses (e.g. crates/catallax-features/tests/saddle_fillets.rs:960-961, 1653-1654) [read].
- Measuring stick, as a library: catallax-compare checks every B-rep display mesh vertex against certified field distance bounds and projects field boundary samples onto B-rep faces; result is a sampled bound, "not a whole-part Hausdorff certificate" (docs/COMPARE.md:7-40) [read]. On-demand only: decision 26 Oct "Field comparison runs on request only" (ROADMAP.md:259) and "not a universal counts-as-built gate" (docs/TARGET-PARTS.md:34) [read].
- In the app: SolidDraft runs an "independent check" that builds each recipe prefix in both kernels and shows "Boundary model" vs "Distance field" status (SD crates/soliddraft-kernel-catallax/src/checks.rs:4-48, 170-254) [read]. Lattice, gyroid, grade field, field mass and mass stand-in are shipped commands (SD crates/soliddraft-app/src/commands.rs:1966-1999) [read]. Print export has a field route (SD crates/soliddraft-app/src/print_check.rs:15) [read].
- Preview: catallax-app exposes the field's PreviewOptions; preview crate uses wgpu (crates/catallax-preview/Cargo.toml) [read]. Roadmap 4.1 "draws from whichever kernel is faster" is not done (GP:1511).

## 2. What works today, measured
- P1-P5 prefix construction at 20M units (docs/M18-NATIVE-FIELDS.md:112-126) [read]: field builds P1 (710 units; B-rep hits budget at slot-cut-6), P2 (3,299), P3 (3,622), P4 (2,894; B-rep 1,204,390), and refuses P5 at the fillet. These are construction costs, "not display costs or a boundary-agreement result". Byte-identical at 1, 2, 4 workers.
- Lattices and gyroids certify STL at 1/10 mm within 2M cells, 1M vertices, 100M units for named design states (FIELD.md:66-80) [read].
- App: gyroid and lattice STL timeouts at 600 s (GP:28, 2022) [read]. Field mass 3.5 to 7.7 s (GP:2604; CX docs/checkpoints/checkpoint-1.md:286) [read].
- Scan fit: three seeded synthetic controls (plate+bore, cone bore, block+boss+fillets) recover planes, cylinders, cones; scan-to-part max 0.04 mm at 0.01 mm noise (M23-SCAN-FIT.md table) [read]. No real scan yet (ROADMAP 4.4 done-when).
- Tests: 466 `#[test]` in catallax-field of 3,706 workspace-wide; 90 in catallax-compare (my grep count).
- No field rows in the checkpoint-2 scoreboard or census that I found (grep of docs/checkpoints/checkpoint-2.md for field/lattice/gyroid: no hits) [read].

## 3. Roles under tolerance-first

| Role | Useful under tolerant B-rep? | Used in shipping products? | Evidence |
|---|---|---|---|
| (a) Independent oracle in tests (inside/outside, volume, mass) | Yes, value rises: once B-rep answers are tolerant, an exact second evaluation is the only cheap independent check. But coverage is weakest exactly where tolerance bites (curved fillets, splines, imported spline STEP). | Not as a second kernel. Kernel vendors test against reference data and other kernels [unread/memory]. | CX features tests (above); GP:1762 "value grows ... check a fitted numerical face"; ARCHITECTURE.md:465 via GP [read] |
| (b) Robust point classification, winding numbers that tolerate gaps | Partly. Field ray parity needs a closed, valid mesh (FIELD.md:11), so it does not tolerate gaps. Generalised winding numbers (Jacobson et al. 2013; fast version Barill et al. 2018) do, and are a separate algorithm, not an SDF. | Yes for GWN in libigl and in mesh tools [unread/memory]. OCCT and Parasolid classify with tolerances on the B-rep itself [unread/memory]. | FIELD.md:11 [read] |
| (c) Lattices, gyroids, TPMS, graded density, many-body blends, complex shells | Yes, unchanged. This is the one role no tolerant B-rep covers well. | Yes: nTop (implicit core), Fusion volumetric lattice plus volumetric API that defines SDFs, Parasolid v38 implicit pattern lattices and v39 "B-reps, facets and implicit geometry within the same model". | Parasolid v39 blog [snippet] https://blogs.sw.siemens.com/plm-components/parasolid-v39-0-advancing-the-future-of-geometry-modeling/ ; Fusion help [snippet] https://help.autodesk.com/cloudhelp/ENU/Fusion-Model/files/SLD-CREATE-VOLUMETRIC-LATTICE.htm ; https://3dprintingindustry.com/news/autodesk-previews-volumetric-modeling-api-for-fusion-240899/ [snippet] |
| (d) 3D printing export (certified STL/3MF) | Yes for field-only parts. For ordinary parts the B-rep mesh route is simpler (printing::brep already exists). | Yes: nTop, Fusion additive, Materialise [unread/memory]. A certified distance bound in the file is unusual. | FIELD.md:23-39 [read]; GP:1764 |
| (e) Collision, clearance, interference | Yes, cheap and robust for assemblies; PRD assigns it to the field (GP:3797, 3911). | Yes in physics and robotics (SDF collision in PhysX 5, Isaac) [unread/memory]; MCAD interference is usually B-rep [unread/memory]. | M18-FINITE-CLEARANCE.md exists [read heading]; field assembly search PR 506 (GP:1510) |
| (f) Healing or rebuilding bad imports (voxel/SDF remesh) | Plausible, not built. Current mesh import refuses invalid meshes. | Yes: OpenVDB voxel remesh in Blender, shrink-wrap in Magics, nTop mesh-to-implicit [unread/memory]. | FIELD.md:11 [read] |
| (g) Preview and GPU display | Weak. B-rep already gives fast unproved previews (catallax-app preview_body uses catallax_mesh). Field preview only wins for field-only features. | Lumen mesh distance fields in games [unread/memory]; not MCAD display. | crates/catallax-app/src/preview.rs:1-60 [read] |
| (h) Simulation (immersed methods) | Yes as a later product role; needs only robust inside/outside and distance, which the field has. Not implemented ("FEA ... not implemented", FIELD.md:7). | Yes: Intact.Simulation (Immersed Method of Moments) on nTop implicits; SimScale immersed-boundary CHT on nTop implicits; nTop Fluids LBM. | https://intact-solutions.com/simulation-built-for-implicit-design/ [snippet]; https://www.simscale.com/press/simscale-and-ntop-integration-announcement/ [snippet] |

Csanady (founder of Shapr3D): implicits are "the strongest known answer to B-rep's math problems and the weakest answer to its workflow problems"; easy to make gyroids and lattices, hard to make an iPhone (csanady-geometry-kernels-summary.md:37-39, 58) [read]. Matches the evidence: the field does well on lattices and prismatic recipes, poorly on the curved fillets that make real parts.

## 4. Cost
- Code: catallax-field 76,643 lines (53,591 src, ~20,569 tests), second-largest crate of 602,736 Rust lines (12.7%). Plus catallax-compare 6,069 and catallax-preview 3,825 (my wc counts).
- Activity: 134 of 793 commits since 20 Oct touch catallax-field (17%); 48 of 451 since 10 Nov (11%) (git log counts).
- Roadmap: all of phase 4 (4.1 to 4.4; 4.5 parked) (ROADMAP.md:189-201). Roadmap says phase 4 "waits until the gaps are closed" (GP:1168), yet all four items are "started" (GP:1511).
- Coupling: catallax-features uses the field only as a dev-dependency, so it does not block B-rep builds. Runtime dependents: app (preview, assembly mass), cli, compare, preview. The real coupling is the PRD rule that every recipe operation has a field version (GP:3789, 3895-3904, 3825 "field version of every feature"), and 4.2 "Every part has both descriptions". Under that rule each new B-rep feature (curved fillets, spline offsets, lofts) carries a second implementation; the P5 verdict shows a single fillet case needs new envelope regions in the field (M18-P5-VERDICT.md:15-21).
- The field and its oracle tests are strongest on analytic recipes, which are the cases the B-rep already handles exactly. Under a tolerant model the risky cases (curved blends, imported splines) are where the field refuses.

## 5. Options (inference)
1. Keep as is (dual-kernel parity, 4.2 "field catches up"). Pros: keeps the PRD promise and all roles open. Cons: every B-rep feature costs twice; the parity goal chases the hardest B-rep work (curved fillets) with no user payoff; competes with phase 1-3 gaps.
2. Narrow to named roles (recommended): (i) test oracle for analytic recipes, frozen at current coverage, with refusals accepted; (ii) field-only features: lattice, gyroid, grades, their certified print and mass, stand-in STEP; (iii) clearance/interference queries. Drop the "every operation in both kernels" rule and 4.2 parity. Field-only features attach to B-rep regions as terminal steps, as they already do (stage refuses CAD after the field, FIELD.md:74). Defer immersed FEA and import healing until a customer needs them. Fix the 600 s STL timeouts because that is a shipped feature failing. Pros: keeps what is unique and shipped, stops the double-implementation tax. Cons: loses the full-part independent check for curved features; "booleans that never fail" fallback goes away.
3. Park (freeze code, no new work). Pros: zero roadmap cost. Cons: shipped lattice/gyroid commands rot; 600 s timeouts stay.
4. Delete. Pros: removes 13% of code. Cons: loses the only differentiator that Parasolid itself is now adding (implicit lattices, v38/v39), loses the test oracle and SolidDraft's lattice/gyroid features.

Recommendation (inference): option 2. Parasolid and Fusion are moving toward mixed B-rep + implicit models for lattices and AM, so "no SDF" is not true of Parasolid in 2025-2026; the industry pattern is B-rep master plus implicit regions, not two full kernels.
