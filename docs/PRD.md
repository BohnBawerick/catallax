# Catallax kernel PRD

Published 22 October 2026. The [research notes](research/) give the sources behind this direction.

[ROADMAP.md](../ROADMAP.md) is authoritative wherever it conflicts with this PRD, by the owner's decision of 15 November 2026. It sets the order, checkpoints, current scope and recorded status. This PRD retains compatible product requirements; it does not claim that planned work is complete.

The [26 October 2026 decision](ARCHITECTURE.md#preview-and-deferred-proof-decision-of-2026-10-26) separates preview/display from deferred proof and makes second-kernel comparison on demand. This changes the build direction, not the current public API or what counts as proven.

The [26 October 2026 print decision](ARCHITECTURE.md#print-accuracy-and-routes-decision-of-2026-10-26) keeps 0.01 mm as Catallax's default print accuracy, uses proven B-rep tessellation for ordinary parts and adaptive certified field cells for field parts, and requires visible output settings and explicit resource trade-offs. It does not claim those routes and reports are already implemented.

## 1. What Catallax is

Catallax is an open-source geometry kernel under the MIT licence.
It is the kernel for SolidDraft, designed for the features SolidDraft's users drive as it grows toward the most usable of SolidWorks's parametric modelling and nTopology's implicit modelling. Anyone may use it; generalising it for others comes after it carries one application well.
It serves the mechanical industry: machined parts, printed parts, drawings, STEP exchange with other CAD tools, FEA and large assemblies.

Catallax was published from the existing B-rep implementation without its Git history. Further fixes arrive as reviewed Catallax changes; current scope and status are recorded in ROADMAP.md.

## 2. The model

- One recipe is the master. The recipe is the feature tree: sketches, extrudes, holes, booleans, fillets, lattices and the rest.
- Two kernels evaluate the same recipe definitions when needed. The B-rep kernel builds exact faces and edges. The field kernel builds a signed distance field, owns field-only features and supplies an on-demand second opinion on B-rep parts.
- Each operation has one mathematical definition, and both kernels use it. A fillet, for example, is the surface swept by a ball of radius r that rolls while touching faces A and B. Two kernels can only agree on operations that have a written definition.
- Faces and edges are named by the operation that made them, for example "the edge between the side face of extrude 1 and the floor of cut 2". A name never depends on the order of faces in a result.
- Nothing is converted from one kernel to the other unless a feature or an output needs it.

| Kernel | Jobs |
|---|---|
| B-rep | STEP export, ordinary-part printing, drawings, CAM, exact dimensions, exact mass properties, mates |
| Field | FEA without meshing, field features for printing such as lattices, gyroids and graded density, interference and clearance checks, booleans under the applicable coverage and result contracts |

## 3. Rules

### 3.1 Agreement

The two kernels agree exactly on exact regions. Permitted approximations carry the proved bounds below.
Geometry creation uses one modelling tolerance of 1 nm, which is 0.000001 mm, for sameness and decisions.
For results published as proven, a tolerance is allowed only in the following two cases, and both are cases of necessity. The separate numerical route below does not widen these proven-result exceptions or their caps.

1. Shapes needing a different representation in the B-rep or in STEP, including general curved-face offsets and fillets. Both kernels retain the exact definition. The B-rep may carry a certified polynomial or rational B-spline surface beside it. Each approximated face records its proved local tolerance: 100 nm target and 1 µm hard cap for kernel-made approximations; 10 µm for imports unless the user explicitly accepts and records a looser limit for that file. A local tolerance is also at most a tenth of the shortest neighbouring edge or narrowest neighbouring face. The body's certificate proves exact agreement except for the listed approximated/imported faces and bounded decisions, with complete coverage, orientation and valid topology. A permitted face approximation does not also have to meet a universal 1 nm whole-result bound. Later operations remeasure against the retained exact definition, so tolerances do not grow by inheritance. Display and export record separate output errors and accuracy requirements. Field-only features remain field-only unless an explicit conversion meets its contract.
2. A decision that takes too long to prove. This is allowed only after the B-rep kernel has proven that the gap behind the decision is smaller than 1 nm. exact-brep already has this outcome, "not exact, with a proven bound". The 1 nm limit is new.

Keep sameness decisions consistent, retain each local tolerance, and check topology before publication as proven. Never clamp a measured tolerance to its cap. [ARCHITECTURE.md](ARCHITECTURE.md#composable-geometric-error-decision-of-2026-10-22) defines composition and the measurement, comparison and STEP hand-offs for these scoped certificates. Section 3.3 permits an explicitly unproved preview before that proof completes.

Where exact proof is unavailable, [ROADMAP milestone 1.1](../ROADMAP.md#plan) permits separately marked numerical geometry under a stated tolerance fixed before the operation. Retain the generating definition, completed checks with their actual scope and sample locations, and any measured repair history. Numerical results remain separate from proven results and from bounded answers that require a complete bound over their declared scope. Zero sampled residuals do not promote a numerical result to either category. Local pointwise distance enclosures, sampled maxima and parameter-space trim checks do not establish a whole-face or whole-body bound, global self-intersection freedom or separation. Report only the evidence actually established; do not add repair movements into an asserted bound on the current body.

The landed [checked numerical build contract](interface-changes/2026-11-15-checked-numerical-build-results.md) and [trim/sewing extension](interface-changes/2026-11-14-numerical-trim-sewing.md) provide separate numerical bodies and evidence within their documented families. Their producers report `exact: false`; [producer guarantees](interface-changes/2026-11-15-numerical-producer-guarantees.md) do not prevent every explicit misuse of generic answer constructors. They grant no exact checkpoint or existing STEP-export authority. Numerical STEP output is a roadmap requirement: an implemented numerical writer must validate its output and mark it approximate before delivery. N01/N02 do not implement that writer; this requirement does not qualify N03 or any broader approximate export. Numerical admission does not change the existing modelling tolerances, proven local caps or work allowances.

Display and export accuracy are output settings, not modelling tolerances. Print accuracy defaults to 0.01 mm and may be changed explicitly for an output under section 5.1. A chosen coarser output does not relax the model's 1 nm decisions, approximation target or local caps.

### 3.2 Refusals

The kernel refuses only geometry that is impossible, the way SolidWorks refuses zero-thickness geometry.
Examples from exact-brep today are `boolean.self-touching-shell`, `shapes.non-positive-radius` and `shapes.torus-self-intersecting`.

Every other refusal means the kernel is not finished.
Examples from exact-brep today:

| Code | What triggers it |
|---|---|
| `intersect.traced-pair-required` | A general intersecting cone or torus pose outside the current closed-form paths |
| `boolean.imprint-ambiguous` | Two pipes of the same diameter crossing at right angles |
| `numbers.operation-budget` | A request exhausting its fixed deterministic work allowance |
| `measure.accuracy-not-reached` | A volume or centroid that could not be proven to the requested accuracy |

About 20 codes of this second kind exist in the code today, including `boolean.classify-ambiguous`, `boolean.merge-ambiguous`, `numbers.budget-exhausted` and `numbers.precision-limit`.
Until the kernel covers such a case, the user sees a failure like a SolidWorks rebuild error, and the case goes on the list of kernel bugs to fix. A print refusal also explains the limiting stage and the accuracy/resource trade-off under section 5.1; it must not be a bare failure or claim the geometry is impossible merely because a cap was reached.
The user never sees a half-built state presented as a proven result and never needs to know a kernel exists. If proof of an already visible preview fails, the app marks it failed and shows the typed refusal for that edit.

Report refusal latency against its applicable preview or background-proof target. Missing a local tolerance cap is unfinished coverage, not impossible geometry. A refusal that misses its clock target is also a performance bug; record it while clock acceptance waits for Phase 5. Fix hangs now rather than deferring them with speed acceptance.

### 3.3 Display

Both kernels act on the same geometry, imported parts included. The faster kernel does the job within its supported scope and result contract; the B-rep owns STEP export and drawings, and the field owns lattices, gyroids and graded density. A requested output or on-demand comparison may require the other representation. Unsupported routes still report their named refusals; selecting a faster route does not promote numerical results or previews to proven geometry.

Preview/display has its own budget for construction preview, display mesh, edge polylines and names needed for selection. The targets are 100 ms for stock solids, 500 ms for booleans and 1,000 ms for a recipe-feature click, measured through the visible selectable edit. A separate background proof budget covers certificates, whole-face deviation bounds, exact requested measurements, export validation and requested field comparison. Its target remains 2,000 ms for the complete requested background job. Proof never blocks the preview click and may be postponed. These click- and proof-clock targets become acceptance gates in [Phase 5](../ROADMAP.md#plan); earlier capability acceptance retains the correctness and output gates in section 10.

A preview uses the same immutable recipe revision and exact definitions as the eventual proof, but need not use the same representation or construction route. It is clearly labelled `unproved` until proof confirms it or replaces it with the proven result. Selection names belong to that revision and are resolved against the completed result; they are not proof. A newer edit invalidates the earlier pairing, and a late proof cannot overwrite it. A previous body may remain visible as the previous revision, not as a preview of new geometry. On model-proof failure, mark the preview `failed` with the typed refusal. Exports, requested measurements and anything published as proven still wait for their applicable proof. An output-specific failure does not invalidate an already proved model. A part built numerically within its stated tolerance shows as finished with its numerical status and completed checks, distinct from proven or bounded status. It follows section 3.1 and its own output validation; a preview alone never grants numerical export authority.

Keep the 1M scoreboard, 20M application and 100M verification-request allowances unchanged. Separating preview and proof must not reset, duplicate or enlarge a caller's deterministic allowance. All modelling tolerances, approximation targets and caps stay unchanged. Display/export accuracy remains an output setting; section 5.1 governs explicit print choices. A proven answer is never wrong.

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

Each approximated face published as proven still needs whole-face deviation bounds against its own exact definition. Numerical evidence follows the separate scope rules in section 3.1. Two-sided per-face distance routines may supply that local proof without rebuilding the whole part in the second kernel. A sampled comparison never substitutes for it. If proof or a requested comparison finds a model failure, the edited feature and its preview are marked failed with the typed refusal, and the case goes on the kernel bug list.

The app never checks by asking the B-rep kernel point by point.
In an earlier measured run exact-brep took 1 to 25 ms per point on ordinary parts, so a useful check of that kind would take 10 seconds to 4 minutes per edit.

### 3.6 Jev

Jev is TypeSafe's decision model. It returns choices and probabilities from text or JSON input.
Catallax uses it only for imports and for choosing intent among options the mathematics has already measured.
It is never used for inside or outside decisions, or for any other decision that has a provable answer.
Every Jev answer is saved in the part, so a rebuild gives the same result and works offline.
Recommended: Jev is an optional helper in SolidDraft, not a dependency of the MIT kernel.

### 3.7 Kept from exact-brep

- Exact results stay identical on every machine and at every thread count. Only explicitly unproved previews that never export may vary in their last digit; this exception grants no variation in proven answers or authoritative geometry.
- No unsafe code.
- Derived geometry is stored by its definition, and coordinates are only a cache.
- A proven answer is never wrong.

### 3.8 Target

Reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels.
The former 1-in-10,000 programme and its grid-cell product goals are retired. Target parts that rebuild in SolidDraft with no wrong answer are milestones toward the target and the published scoreboard. Click-speed and proof-clock acceptance wait for Phase 5; correctness, completed proof, export, naming and application gates remain in force.

### 3.9 Threads

Rust's safe threads run rebuilds and checks in the background, so the user is never blocked. Ordinary hardware or GPU arithmetic may vary in the last digit only for the explicitly unproved, non-exportable previews described in section 3.7; exact results retain full determinism.

## 4. Features

This table describes product roles, not completion. Supported families and remaining work follow ROADMAP's “Where we are” and milestones; exact, numerical, bounded and preview statuses retain section 3.1's distinct claims.

| Feature | B-rep kernel | Field kernel |
|---|---|---|
| Sketch, extrude, revolve | Exact | Exact, from a formula |
| Hole, countersink | Exact where supported; remaining general holes follow milestone 1.4 | Exact |
| Boolean | Exact. This is exact-brep's core | Exact where supported; unsupported coverage remains a named refusal |
| Fillet or chamfer on an edge between faces the recipe defines | Exact for flat to flat and for a flat face meeting a perpendicular hole. Curved fillets retain their definition and may use a fitted B-rep carrier under section 3.1; general mixed and spline families follow milestones 1.4 and 1.5 | Exact rolling ball, because the recipe names both faces. Slower between two curved faces |
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
| Display of one part | Whichever kernel is faster within its supported scope and status contract | No |
| Display of a big assembly | Meshes with level of detail, from either kernel | No |
| Picking faces and edges | B-rep faces, or recipe-defined faces in the field | No. Freeform field surfaces have nothing to pick |
| STEP | B-rep, at the accuracy chosen at export | Yes, when the part has field features |
| Drawing | B-rep | Yes, when the part is printed from the field |
| CAM | B-rep, through STEP | No |
| FEA | Field | Yes. It must match the exported part |
| Print | Ordinary parts from proven B-rep faces; field parts from adaptive certified field cells, at the chosen print accuracy | Yes when a field print is inspected against a B-rep drawing or the output otherwise combines representations; no automatic field comparison for B-rep-only printing |
| Mass, centre of gravity | B-rep exact, plus the field for field-only features | Yes. The user sees one number |
| Mates | B-rep planes and axes, or recipe-defined field faces | Yes, together with interference |
| Interference, clearance | Field | Yes, together with mates |

### 5.1 Print accuracy and resource choices

Default to 0.01 mm print accuracy and a visible 2,000,000-cell cap. Keep the separate 1,000,000-vertex default. These are output controls, recorded even when supplied by default. Report the cell cap on every successful or refused print, with its applicability and actual usage; a B-rep route without cell traversal does not consume field cells.

Ordinary parts print by tessellating proven B-rep faces and trims within the chosen deviation. Flat polygons need no accuracy-driven interior refinement; curves need enough subdivision to meet their bound. Prove file-coordinate rounding and composition of model, tessellation and encoding error, manifold vertex links and cavity orientation. A display mesh, preview or successful file write is not automatically a print.

Lattices, gyroids, graded density and other field parts use adaptive cells, fine near detail and coarse on flat or gentle regions only when their local certificate permits it. Keep certified error per cell size, complete two-sided boundary coverage and valid oriented topology across refinement levels. Neither route may delete small components or cavities, silently switch representations or publish a partial mesh. A B-rep print does not satisfy a native-field acceptance case.

When the request cannot fit, name the limiting stage and report requested accuracy, effective caps, actual cells/vertices/facets and work used. Give labelled estimates of cells, triangles, peak memory and time at that accuracy, a coarser accuracy estimated to fit the current limits, and the higher-cap cost of retaining accuracy where supported. Explain when no fitting candidate or credible cost estimate is available. Identify any remaining vertex or work limit. Bounded coarse traversal, live-storage accounting and identified throughput calibration underpin these estimates; [the architecture decision](ARCHITECTURE.md#visible-limits-and-actionable-refusals) specifies their scope and uncertainty. Estimates are not proof or a promise of completion.

The user explicitly chooses a new request; controls never weaken silently. A chosen 0.05 mm print must prove 0.05 mm and must not be reported as a 0.01 mm print. Successful files carry the chosen accuracy and cap in documented metadata, with a byte-bound manifest retaining the requested settings, achieved file-to-definition bound, units, encoding/rounding, preservation policy and completeness/topology/orientation status. For STL this needs an explicit header/consumer contract because STL standardizes neither units nor certified accuracy. A format that cannot retain the settings losslessly must not silently discard them.

The 1 nm modelling tolerance, 100 nm approximation target, 1 µm kernel-made cap and other local limits remain unchanged, as do the distinct 1M/20M/100M work allowances. Raising a print cell cap does not raise a work allowance. Existing native-field qualification controls remain in force.

## 6. What the kernel provides

1. A recipe language with one mathematical definition per operation, and names given by the operation that made each face and edge.
2. A B-rep version of every B-rep-native feature, with the necessity surfaces from section 3.1 held exactly by their definition.
3. A field version of every feature, including imports, with the recipe's faces kept selectable.
4. Conversion in both directions on demand, only when it meets its contract. A field-only feature stays field-only if conversion cannot meet that contract. Field to B-rep conversion at a chosen accuracy and the section 7 stand-in are distinct outputs; a stand-in is not full conversion of the field geometry. B-rep to field follows the same scoped hand-off requirements.
5. Meshes from both kernels at a chosen accuracy.
6. The on-demand surface comparison and requested hand-off checks from section 3.5.
7. The queries FEA and interference checks need: inside or outside, distance and surface direction. Exact mass from the B-rep side.
8. Background rebuilds and checks on several threads.

## 7. Export and import

- **STEP accuracy.** The user chooses the accuracy level at export.
- **Stand-in for very large parts.** When a full conversion would be huge, STEP and drawings use a stand-in the way nTop does it: the outer envelope plus bodies that match the mass, centre of gravity and inertia. nTop's example lattice cube was over 600 MB as STEP and under 30 kB in its own format. nTop suggests full conversion only up to about 10,000 beams or 1,000 gyroid cells. The size that switches to the stand-in is not decided yet.
- **Printing.** Slices, 3MF or STL under section 5.1. Ordinary parts use a print-qualified tessellation of the proven B-rep; field parts use adaptive certified cells. The chosen accuracy and cap travel with the file, and an output that cannot fit returns an explicit trade-off for the user to choose.
- **STEP import.** The imported body enters the recipe as a fixed piece. The field kernel computes its field from it. Jev decides whether small gaps and offsets are intended or noise from the system that wrote the file.
- **Mesh and scan import.** Aim for a clean, workable part with the smallest deviation, as B-rep and as field, without requiring a reconstructed recipe or feature tree. Use classic fitting first; use learned segmentation such as DualBrep only if classic fitting is insufficient. Never silently treat a scan as a solid. Report deviation in both directions, its measurement scope and whether it is sampled evidence or a complete bound under section 3.1. Matching RECIPE section 4.19.1 wording requires a separate coordinated change; this decision does not perform a broad recipe migration.

## 8. Outside this PRD for now

- **FEA.** Remains parked under [ROADMAP milestone 4.5](../ROADMAP.md#phase-4-field-lattices-printing-and-scans) until the owner restarts it. The field kernel gives the solver inside or outside, distance and surface direction. How the solver works gets its own discussion. The historical FEA-first order is not adopted.
- Rotation and rendering throughput beyond the preview/display targets in section 3.3. Measure those separately before deciding whether more work is needed.
- **Large-assembly display.** Research found that CAD tools draw big assemblies as meshes with level of detail, lazy loading and instancing. No CAD viewer was found that draws assemblies from SDFs. In assemblies the field kernel's role is interference and clearance.

## 9. How this relates to existing work

### 9.1 exact-brep becomes the B-rep half

At publication, the [Catallax repository](https://github.com/BohnBawerick/catallax) received the B-rep code without history, with renamed workspace and preserved exchange formats. Further source fixes arrive as reviewed Catallax pull requests, not an implicit sync. The field kernel and recipe language were new work at that point; ROADMAP records their current scope.

The original publication carried over the following from exact-brep, about 32,000 lines of Rust. This is historical provenance, not a current capability inventory:

- exact numbers and certified decisions with work budgets;
- plane, cylinder, cone, sphere and torus builders, with shapes stored by definition;
- booleans on plane, cylinder and sphere pairs, through the five stages intersect, imprint, classify, keep and sew;
- proven volume, area and centroid;
- STEP AP242 export;
- the same answer on every machine and thread count, and no unsafe code;
- the command-line interface the judge uses.

What changes in exact-brep's design documents (`docs/DESIGN.md`, `docs/ARCHITECTURE.md`) and in [PLAN.md](archive/PLAN-2026-11.md):

- what it is part of: the B-rep half of Catallax;
- the refusal rule in DESIGN R1, per section 3.2;
- the proven-bound decision outcome in R1 and R3, with geometric gaps below 1 nm per section 3.1;
- the target, per section 3.8, which retires the 1-in-10,000 programme and its grid cells as the goal;
- speed per click now matters, because the kernel serves an interactive CAD tool;
- faces and edges named by the operation that made them.

The original publication listed these gaps. Several now have supported families in ROADMAP's “Where we are”; this historical list must not be read as current absence or completion of the remaining general cases:

- general cones and tori in booleans beyond the closed-form and coaxial paths;
- sketch-based features: profile, extrude, revolve;
- fillets and chamfers;
- shell, offset and draft;
- sweeps and lofts;
- patterns and mirror;
- meshes for display;
- STEP import;
- seam curves in STEP for interoperability with other CAD tools. The historical exchange forbade them; current consumers must follow their applicable versioned seam contract;
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

SolidDraft uses Catallax only; its kernel slot is not swappable. Integrate each kernel operation into the app as it lands, using minimal adapters for the shared operation definitions, persistent names and hand-off checks. No broad recipe-format migration is authorised by this decision. Coordinate any necessary RECIPE section 6 wording separately. The kernel remains usable independently of SolidDraft.

## 10. Direction, in order

Follow [ROADMAP's phase order and checkpoints](../ROADMAP.md#plan):

1. Phase 0 closes the release and its evidence without claiming that preparation alone qualifies it.
2. Phase 1 closes numerical modelling under stated tolerances, real-file import, spline cuts and joins, everyday edge/wall families and spline fillets. Integrate each operation into SolidDraft as it lands, browser first, with section 9.4's minimal adapters and unchanged result gates.
3. Phase 2 adds the NURBS toolkit, surface bodies, general sweeps and direct edits in the roadmap's order and dependencies.
4. Phase 3 covers drawing views, assemblies/exchange, clash/clearance/mass and the remaining modelling lines.
5. Phase 4 covers field preview, field parity, lattices/printing and imports/scans in both kernels. FEA remains parked until the owner restarts it at a checkpoint. Field parity never blocks a narrow B-rep feature; requested comparisons and hand-offs retain section 3.5's coverage requirements.
6. Phase 5 addresses speed and simplification, GPU work, Windows packages and signed downloads. Use every core by design with deterministic result assembly; clock acceptance waits for this phase. App-alongside and standing lanes follow ROADMAP throughout.

At checkpoints, stop new work, merge or park running work with its limits recorded, gather the scoreboard/census/refusals and update the roadmap; the owner chooses the next phase's order. A wrong proven answer stops work on that method immediately.

Target-part and milestone acceptance follows the [current roadmap](../ROADMAP.md#plan) and [target-part contracts](TARGET-PARTS.md). Capability milestones require a completed build within the applicable fixed work allowance, with its declared result category and evidence, and the existing correctness, proof, export, name-survival and application checks. Numerical status does not satisfy a gate requiring a proven result, a whole-body bound or a validated export. Click- and proof-clock acceptance wait for Phase 5; proof completion itself does not wait. List field comparison as an on-demand check with its coverage and outcome, not a universal landing or built gate. Preview alone does not count as a built part.

Cheaper-proof and simplification work belongs in Phase 5, in the measured order in ROADMAP, retaining proofs from construction, prepared bounds, reuse and coarse-to-fine refinement where valid. No independent body checker, proof verifier, durable proof-job machinery or bounded-modelling contract project is scheduled: ROADMAP closes or rejects those projects. Existing proof, output-validation and evidence-scope requirements remain obligations of any result that claims them; this does not reopen a separate certificate-system build.

## 11. Decisions and remaining open items

The obsolete 50-second individual clock is not an acceptance rule. The existing fixed deterministic work allowance selects verification outcomes, not machine load or elapsed time. The application edit allowance stays 20M; retain the separate historical 1M scoreboard and 100M verification-request contracts. Clock time never fails a change before Phase 5; record latency and fix hangs now. This documentation decision does not claim that an implementation has changed.

The tolerance decision of 23 October (proposal 3) and on-demand comparison decision of 26 October (proposal 5) were already applied in sections 3.1 and 3.5; they are not new decisions or duplicate implementation work. The separately marked numerical route above is the subsequent scoped direction.

Remaining items:

- The size at which export switches to the stand-in.
- Remaining recipe-language and rebuild-stable-name requirements follow the applicable roadmap milestones; the recipe language itself already exists.
- Whether the field kernel builds on Fidget, a Rust field engine under MPL-2.0, or is written from scratch. This was not discussed. MPL-2.0 keeps its own files under MPL inside an MIT project.
- The FEA design.
- The on-screen speed measurement.
- Any unsupported exact/exception evidence needed by a scoped result or output must be reported honestly under section 3.1. This is not an open project to build the certificate systems closed by ROADMAP.

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
