# Catallax kernel PRD

Published 22 October 2026. The [research notes](research/) give the sources behind this direction.

This document sets the direction for Catallax. It describes requirements and planned features, not a claim that they are already implemented.

The [26 October 2026 decision](ARCHITECTURE.md#preview-and-deferred-proof-decision-of-2026-10-26) separates preview/display from deferred proof and makes second-kernel comparison on demand. This changes the build direction, not the current public API or what counts as proven.

## 1. What Catallax is

Catallax is an open-source geometry kernel under the MIT licence.
It is the kernel for SolidDraft and for anyone else who wants it.
It serves the mechanical industry: machined parts, printed parts, drawings, STEP exchange with other CAD tools, FEA and large assemblies.

Catallax is a new public repository, exported from the existing B-rep implementation without its Git history. The original repository remains the source for separately reviewed fixes during the transition.

## 2. The model

- One recipe is the master. The recipe is the feature tree: sketches, extrudes, holes, booleans, fillets, lattices and the rest.
- Two kernels evaluate the same recipe definitions when needed. The B-rep kernel builds exact faces and edges. The field kernel builds a signed distance field, owns field-only features and supplies an on-demand second opinion on B-rep parts.
- Each operation has one mathematical definition, and both kernels use it. A fillet, for example, is the surface swept by a ball of radius r that rolls while touching faces A and B. Two kernels can only agree on operations that have a written definition.
- Faces and edges are named by the operation that made them, for example "the edge between the side face of extrude 1 and the floor of cut 2". A name never depends on the order of faces in a result.
- Nothing is converted from one kernel to the other unless a feature or an output needs it.

| Kernel | Jobs |
|---|---|
| B-rep | STEP export, drawings, CAM, exact dimensions, exact mass properties, mates |
| Field | FEA without meshing, field features for printing such as lattices, gyroids and graded density, interference and clearance checks, booleans that never fail |

## 3. Rules

### 3.1 Agreement

The two kernels agree exactly on exact regions. Permitted approximations carry the proved bounds below.
Geometry creation uses one modelling tolerance of 1 nm, which is 0.000001 mm, for sameness and decisions.
A tolerance is allowed only in two cases, and both are cases of necessity.

1. Shapes needing a different representation in the B-rep or in STEP, including general curved-face offsets and fillets. Both kernels retain the exact definition. The B-rep may carry a certified polynomial or rational B-spline surface beside it. Each approximated face records its proved local tolerance: 100 nm target and 1 µm hard cap for kernel-made approximations; 10 µm for imports unless the user explicitly accepts and records a looser limit for that file. A local tolerance is also at most a tenth of the shortest neighbouring edge or narrowest neighbouring face. The body's certificate proves exact agreement except for the listed approximated/imported faces and bounded decisions, with complete coverage, orientation and valid topology. A permitted face approximation does not also have to meet a universal 1 nm whole-result bound. Later operations remeasure against the retained exact definition, so tolerances do not grow by inheritance. Display and export record separate output errors and accuracy requirements. Field-only features remain field-only unless an explicit conversion meets its contract.
2. A decision that takes too long to prove. This is allowed only after the B-rep kernel has proven that the gap behind the decision is smaller than 1 nm. exact-brep already has this outcome, "not exact, with a proven bound". The 1 nm limit is new.

Decide once and reuse the record; merge objects declared the same; retain each local tolerance; and check topology before publication as proven. Never clamp a measured tolerance to its cap. [ARCHITECTURE.md](ARCHITECTURE.md#composable-geometric-error-decision-of-2026-10-22) defines composition and the measurement, comparison and STEP hand-offs for these scoped certificates. Section 3.3 permits an explicitly unproved preview before that proof completes.

### 3.2 Refusals

The kernel refuses only geometry that is impossible, the way SolidWorks refuses zero-thickness geometry.
Examples from exact-brep today are `boolean.self-touching-shell`, `shapes.non-positive-radius` and `shapes.torus-self-intersecting`.

Every other refusal means the kernel is not finished.
Examples from exact-brep today:

| Code | What triggers it |
|---|---|
| `intersect.traced-pair-required` | A general intersecting cone or torus pose outside the current closed-form paths |
| `boolean.imprint-ambiguous` | Two pipes of the same diameter crossing at right angles |
| `numbers.deadline` | A whole operation running past 50 seconds, such as the near-miss cases |
| `measure.accuracy-not-reached` | A volume or centroid that could not be proven to the requested accuracy |

About 20 codes of this second kind exist in the code today, including `boolean.classify-ambiguous`, `boolean.merge-ambiguous`, `numbers.budget-exhausted` and `numbers.precision-limit`.
Until the kernel covers such a case, the user sees a plain failure like a SolidWorks rebuild error, and the case goes on the list of kernel bugs to fix.
The user never sees a half-built state presented as a proven result and never needs to know a kernel exists. If proof of an already visible preview fails, the app marks it failed and shows the typed refusal for that edit.

### 3.3 Display

The kernel that owns a feature builds it and shows it.
Extrudes, holes and fillets belong to the B-rep kernel.
Lattices, gyroids and graded density belong to the field kernel.
The other kernel builds its version only when a requested output or an on-demand comparison needs it.

Preview/display has its own budget for construction preview, display mesh, edge polylines and names needed for selection. The targets are 100 ms for stock solids, 500 ms for booleans and 1,000 ms for a recipe-feature click, measured through the visible selectable edit. A separate background proof budget covers certificates, whole-face deviation bounds, exact requested measurements, export validation and requested field comparison. Its target remains 2,000 ms for the complete requested background job. Proof never blocks the preview click and may be postponed.

A preview uses the same immutable recipe revision and exact definitions as the eventual proof, but need not use the same representation or construction route. It is clearly labelled `unproved` until proof confirms it or replaces it with the proven result. Selection names belong to that revision and are resolved against the completed result; they are not proof. A newer edit invalidates the earlier pairing, and a late proof cannot overwrite it. A previous body may remain visible as the previous revision, not as a preview of new geometry. On model-proof failure, mark the preview `failed` with the typed refusal. Exports, requested measurements and anything published as proven still wait for their applicable proof. An output-specific failure does not invalidate an already proved model.

Keep the 1M scoreboard, 20M application and 100M verification-request allowances unchanged. Separating preview and proof must not reset, duplicate or enlarge a caller's deterministic allowance. All modelling tolerances, approximation targets and caps, and display/export accuracy settings stay unchanged. A proven answer is never wrong.

### 3.4 Hand-off points

Agreement must hold wherever a result from one kernel is used together with the other kernel's version of the same part.

| Hand-off | Why it must agree |
|---|---|
| FEA on the field, then STEP export from the B-rep | The stress result must describe the part the shop machines |
| Mates on B-rep faces, together with interference checks in the field | The field part must sit where the mate put it |
| A drawing from the B-rep of a part printed from the field | The dimensions on paper must describe the printed part |
| A B-rep feature that meets a field feature, such as a fillet where a plane meets a gyroid, or a lattice inside a machined shell | Where the two meet, both kernels need the same boundary |
| Mass and centre of gravity | The user sees one number |

Agreement does not matter for a picture on screen, or for a workflow that only ever uses one kernel, such as a machined part going straight to STEP or a lattice going straight to a printer.

### 3.5 Checks

The field's second-kernel comparison runs only when the user or a verification run asks for it. It does not run on every edit and is not a prerequisite for a B-rep feature or its landing.

1. In development, explicit verification runs use the independent tester, exact and certified checks using the Raghothama and Shapiro method, and random recipes through both kernels. Claims of proof need exact agreement or complete coverage of permitted local bounds, including the 1 nm decision limit.
2. In the app, a requested surface comparison runs in the background, separately from preview delivery. It can evaluate the field at B-rep mesh vertices and check the reverse direction. The existing approximately 0.01 mm mesh accuracy makes this a sampled diagnostic, not a whole-face or whole-result proof. Its cost must be measured, not assumed from another engine.
3. Exports and requested measurements wait for their own validation. A user request combining both representations also requests the scoped hand-off agreement proof from section 3.4 before delivery. A single-kernel export does not automatically request a whole-part comparison.

For a requested comparison, the field builds supported steps independently. Where it cannot yet build a step, it may take over the certified B-rep result, explicitly labelling that step non-independent and reporting which later steps it checked independently. Unsupported takeover or missing proof remains unresolved or a typed refusal. The report must not call partial independent coverage full agreement.

Each approximated face still needs whole-face deviation bounds against its own exact definition. Two-sided per-face distance routines may supply that local proof without rebuilding the whole part in the second kernel. A sampled comparison never substitutes for it. If proof or a requested comparison finds a model failure, the edited feature and its preview are marked failed with the typed refusal, and the case goes on the kernel bug list.

The app never checks by asking the B-rep kernel point by point.
In an earlier measured run exact-brep took 1 to 25 ms per point on ordinary parts, so a useful check of that kind would take 10 seconds to 4 minutes per edit.

### 3.6 Jev

Jev is TypeSafe's decision model. It returns choices and probabilities from text or JSON input.
Catallax uses it only for imports and for choosing intent among options the mathematics has already measured.
It is never used for inside or outside decisions, or for any other decision that has a provable answer.
Every Jev answer is saved in the part, so a rebuild gives the same result and works offline.
Recommended: Jev is an optional helper in SolidDraft, not a dependency of the MIT kernel.

### 3.7 Kept from exact-brep

- The same question gives the same answer on every machine and at every thread count.
- No unsafe code.
- Derived geometry is stored by its definition, and coordinates are only a cache.
- A proven answer is never wrong.

### 3.8 Target

Reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels.
The former 1-in-10,000 programme and its grid-cell product goals are retired.

### 3.9 Threads

Rust's safe threads run rebuilds and checks in the background, so the user is never blocked.

## 4. Features

| Feature | B-rep kernel | Field kernel |
|---|---|---|
| Sketch, extrude, revolve | Exact | Exact, from a formula |
| Hole, countersink | Exact. Coaxial cylindrical countersinks now work through booleans; general poses and the named hole feature remain unfinished | Exact |
| Boolean | Exact. This is exact-brep's core | Exact, never fails |
| Fillet or chamfer on an edge between faces the recipe defines | Exact for flat to flat and for a flat face meeting a perpendicular hole. A fillet between two curved faces is held exactly by its definition and becomes a spline only in the STEP file | Exact rolling ball, because the recipe names both faces. Slower between two curved faces |
| Shell, offset, thicken | Exact for planes, cylinders, spheres and tori | Simple only where the field is a true distance. After booleans it is not, so the field measures distance from the faces themselves |
| Draft | Exact. Planes stay planes and cylinders become cones | Exact |
| Sweep, loft | Exact spline faces | Distance to the faces is computed numerically, so it is slower |
| Pattern, mirror | Each copy costs a boolean. One SOLIDWORKS user projected about 3 hours to rebuild a pattern of 1,320 copies | Copies are nearly free |
| Lattice, gyroid, graded density, field-driven thickness | Built only at export, at the chosen accuracy, or as the stand-in described in section 7 | Native |
| Smooth blend | Fitted splines only. There are no edges to pick | Native |
| STEP import | Native | Computed from the faces. Exact formulas for analytic faces, numeric for splines |
| Mesh or scan import | Converted, as described in section 7 | Native, from the mesh |

## 5. Outputs

| Output | Served by | Hand-off point |
|---|---|---|
| Display of one part | The kernel that owns the feature | No |
| Display of a big assembly | Meshes with level of detail, from either kernel | No |
| Picking faces and edges | B-rep faces, or recipe-defined faces in the field | No. Freeform field surfaces have nothing to pick |
| STEP | B-rep, at the accuracy chosen at export | Yes, when the part has field features |
| Drawing | B-rep | Yes, when the part is printed from the field |
| CAM | B-rep, through STEP | No |
| FEA | Field | Yes. It must match the exported part |
| Print | Field slices or 3MF, or the B-rep's mesh | Yes, when the part is inspected against a B-rep drawing |
| Mass, centre of gravity | B-rep exact, plus the field for field-only features | Yes. The user sees one number |
| Mates | B-rep planes and axes, or recipe-defined field faces | Yes, together with interference |
| Interference, clearance | Field | Yes, together with mates |

## 6. What the kernel provides

1. A recipe language with one mathematical definition per operation, and names given by the operation that made each face and edge.
2. A B-rep version of every B-rep-native feature, with the necessity surfaces from section 3.1 held exactly by their definition.
3. A field version of every feature, including imports, with the recipe's faces kept selectable.
4. Conversion in both directions on demand. Field to B-rep at a chosen accuracy, plus the stand-in. B-rep to field.
5. Meshes from both kernels at a chosen accuracy.
6. The on-demand surface comparison and requested hand-off checks from section 3.5.
7. The queries FEA and interference checks need: inside or outside, distance and surface direction. Exact mass from the B-rep side.
8. Background rebuilds and checks on several threads.

## 7. Export and import

- **STEP accuracy.** The user chooses the accuracy level at export.
- **Stand-in for very large parts.** When a full conversion would be huge, STEP and drawings use a stand-in the way nTop does it: the outer envelope plus bodies that match the mass, centre of gravity and inertia. nTop's example lattice cube was over 600 MB as STEP and under 30 kB in its own format. nTop suggests full conversion only up to about 10,000 beams or 1,000 gyroid cells. The size that switches to the stand-in is not decided yet.
- **Printing.** Slices, 3MF or STL, either from the field directly or from the B-rep's mesh.
- **STEP import.** The imported body enters the recipe as a fixed piece. The field kernel computes its field from it. Jev decides whether small gaps and offsets are intended or noise from the system that wrote the file.
- **Mesh and scan import.** An AI model writes a recipe from the mesh, following the CAD-Recode and CADFit research. A numeric fit refines the dimensions. Jev picks the intent among measured options. The B-rep kernel builds the result exactly, and the app measures and shows the deviation. No published method guarantees a tolerance, so the measured deviation is the guarantee.

## 8. Outside this PRD for now

- **FEA.** The field kernel gives the solver inside or outside, distance and surface direction. How the solver works gets its own discussion.
- Rotation and rendering throughput beyond the preview/display targets in section 3.3. Measure those separately before deciding whether more work is needed.
- **Large-assembly display.** Research found that CAD tools draw big assemblies as meshes with level of detail, lazy loading and instancing. No CAD viewer was found that draws assemblies from SDFs. In assemblies the field kernel's role is interference and clearance.

## 9. How this relates to existing work

### 9.1 exact-brep becomes the B-rep half

Use the new public [Catallax repository](https://github.com/BohnBawerick/catallax). Export the existing B-rep code without history, rename the workspace and keep the legacy exchange contract working. Further source fixes arrive as reviewed Catallax pull requests, not an implicit sync. The field kernel and recipe language are new work.

What carries over unchanged from exact-brep, about 32,000 lines of Rust:

- exact numbers and certified decisions with work budgets;
- plane, cylinder, cone, sphere and torus builders, with shapes stored by definition;
- booleans on plane, cylinder and sphere pairs, through the five stages intersect, imprint, classify, keep and sew;
- proven volume, area and centroid;
- STEP AP242 export;
- the same answer on every machine and thread count, and no unsafe code;
- the command-line interface the judge uses.

What changes in exact-brep's design documents (`docs/DESIGN.md`, `docs/ARCHITECTURE.md`) and in [PLAN.md](PLAN.md):

- what it is part of: the B-rep half of Catallax;
- the refusal rule in DESIGN R1, per section 3.2;
- the proven-bound decision outcome in R1 and R3, with geometric gaps below 1 nm per section 3.1;
- the target, per section 3.8, which retires the 1-in-10,000 programme and its grid cells as the goal;
- speed per click now matters, because the kernel serves an interactive CAD tool;
- faces and edges named by the operation that made them.

What the B-rep half still lacks:

- general cones and tori in booleans beyond the closed-form and coaxial paths;
- sketch-based features: profile, extrude, revolve;
- fillets and chamfers;
- shell, offset and draft;
- sweeps and lofts;
- patterns and mirror;
- meshes for display;
- STEP import;
- seam curves in STEP for interoperability with other CAD tools. The existing exchange contract forbids them, so adding them requires a versioned contract change;
- the export accuracy levels and the stand-in.

### 9.2 New pieces

- the recipe language;
- the field kernel;
- the surface comparison and hand-off checks;
- the conversions in both directions;
- the FEA and interference queries;
- Jev on imports;
- the field-only features.

### 9.3 Independent testing

An independent tester built without seeing the kernel's code tests the B-rep half through [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md). Its internal design is outside this repository.

### 9.4 SolidDraft

SolidDraft is the application this kernel serves. Integration must use the shared operation definitions, persistent names and hand-off checks in this PRD. The kernel remains usable independently of that application.

## 10. Direction, in order

1. Publish the clean Catallax repository and align its design, architecture and plan with this PRD.
2. Continue independent B-rep testing and review subsequent fixes without changing the exchange contract in this publication.
3. Define the recipe language. Both kernels and consuming applications need the same definitions and names.
4. Grow the field kernel toward B-rep feature parity and provide the on-demand background comparison. Required cross-kernel hand-offs wait for their scoped proof, but an ordinary B-rep click does not wait for a second build.
5. Widen the kernels along the feature list in section 4 using the shared definitions. A narrow B-rep feature may land with its own proof before field parity; it need not land in both kernels together. Requested comparison reports independent and taken-over steps as section 3.5 requires. Field-only features remain owned and built by the field.
6. Integrate Catallax into SolidDraft once it covers the operations the application needs.
7. Add the field-only features, the export levels and the stand-in, then FEA and interference. FEA is designed in its own discussion first.

Target-part and milestone acceptance follows [PLAN sections 9 and 10.1](PLAN.md#9-target-parts): preview/display meets the interactive targets, background proof completes within its separate budget, and the existing export, name-survival and application checks pass. List field comparison as an on-demand check with its coverage and outcome, not a universal landing or built gate. Preview alone does not count as a built part.

Prioritise cheaper-proof research: proof from recorded construction facts, certificate caching and reuse across edits, proof scoped to requested exports or measurements, and coarse certificates refined on demand. These are research candidates, not decided designs or weaker accuracy requirements.

## 11. Open items

- The size at which export switches to the stand-in.
- The recipe language, including rebuild-stable names.
- Whether the field kernel builds on Fidget, a Rust field engine under MPL-2.0, or is written from scratch. This was not discussed. MPL-2.0 keeps its own files under MPL inside an MIT project.
- The FEA design.
- The on-screen speed measurement.
- The current 50-second `numbers.deadline` cutoff can answer or refuse the same question depending on machine load. This conflicts with the same-answer rule; resolution needs a separate design and contract change, not a code change in this publication.
- Enforcing and composing complete exact/exception certificates over the inherited mark APIs, with the 1 nm decision limit and recorded local face caps.

## Sources

Research notes, each fact with its link and date:

- [research/research-A-commercial-hybrids.md](research/research-A-commercial-hybrids.md): Parasolid, NX, Altair, nTop, Fusion, Dyndrite and others.
- [research/research-B-pure-sdf.md](research/research-B-pure-sdf.md): PicoGK and LEAP 71, Fidget, libfive and other pure field tools.
- [research/research-C-manufacturing-routes.md](research/research-C-manufacturing-routes.md): printing, CAM, mesh-to-B-rep tools, AI reconstruction, STEP and 3MF standards.
- [research/research-D-academic-hybrids.md](research/research-D-academic-hybrids.md): hybrid representations, conversions, fillet theory, Raghothama and Shapiro.
- [research/research-E-implicit-features.md](research/research-E-implicit-features.md): implicit features in practice and how tools mix them with B-rep.
- [research/research-F-large-assemblies.md](research/research-F-large-assemblies.md): large-assembly display, mates, patterns, interference.
- Pratt and Kim, "Experience in the Exchange of Procedural Shape Models using ISO 10303 (STEP)", [NIST paper](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=904157).

exact-brep sources used: `docs/DESIGN.md` R1 to R3, `docs/JUDGE-INTERFACE.md`, `docs/pieces/05-boolean-engine.md`, and the refusal codes in `crates/*/src`.
