# Catallax architecture

A [recipe](RECIPE.md) is the part's feature tree and source of truth. The boundary representation (B-rep) builds faces, edges and solids. The signed distance field evaluates the same operation definitions, owns field-only features such as lattices and gyroids, and can check the B-rep on request. SolidDraft owns documents, commands and user interaction. See the [roadmap](../ROADMAP.md) for planned work and [capabilities](CAPABILITIES.md) for the limits of the current implementation.

## Where changes belong

| Area | Crates and contracts |
| --- | --- |
| Numbers and proof | `catallax-numbers` owns exact and enclosed arithmetic, work accounting, and [decisions](pieces/01-numbers-and-decisions.md). |
| Geometry and topology | `catallax-shapes` stores [geometry definitions](pieces/02-shapes.md). `catallax-topology` stores and [checks faces, loops and bodies](pieces/03-topology.md). |
| Construction | `catallax-intersect` [finds meetings](pieces/04-intersections.md). `catallax-boolean` [imprints and joins or cuts solids](pieces/05-boolean-engine.md). `catallax-features` evaluates the B-rep recipe. |
| Queries and exchange | `catallax-measure` owns [B-rep measurements](pieces/06-measurement.md). `catallax-import` [admits and repairs STEP geometry](pieces/08-step-import-repair.md). `catallax-step` [writes STEP](pieces/07-step-export.md). |
| Field and comparison | `catallax-field` implements [field queries and output](FIELD.md). `catallax-compare` checks [agreement in the requested scope](COMPARE.md). |
| Display and app | `catallax-mesh` builds display meshes with source-face IDs. `catallax-app` maps those IDs to recipe names and manages revision-bound preview and build results. |
| Recipe and CLI | `catallax-recipe` validates and serializes recipe data without evaluating geometry. `catallax-cli` exposes the [legacy exchange](JUDGE-INTERFACE.md) and [deterministic verification exchange](VERIFY-INTERFACE.md). |

## Geometry by definition

Store the support and construction that define a curve, surface or vertex. Treat rounded coordinates and fitted curves as derived values. An edge or face name comes from the operation and its inputs, not iteration order or a temporary handle. Keep names, trims, orientation and the finite domain when a derived support reduces to an analytic one. A new support needs charts, predicates, intersection, integration and export before a caller can treat it as complete.

Use rational arithmetic when it suffices, and retain algebraic expressions and certified enclosures when a legal rational input produces an algebraic result. Record tangency, coincidence and incidence only when construction or proof establishes them. An approximate cache cannot choose a root or topology. Keep curved topology general: a face need not be planar or have one outer loop, and a closed edge need not have vertices.

## Tolerance and proof

- Exact entities have no stored tolerance. Exact decisions keep their existing predicates. For a proven result, preserve exact regions. A near-coincidence needs a proved geometric gap below 1 nm. Certified replacement geometry must retain the original definition, account for error over the whole claimed scope, and respect its local caps. For that route, the fit target is 100 nm and the kernel-made cap is 1 µm; the import cap is 10 µm unless the user opts in to a recorded exception. Also cap the error at a tenth of the smallest neighbouring feature. These are not the separate numerical admission rules. A display mesh, point residual or isolated root enclosure is not a model certificate.
- If exact construction is unavailable or too costly, build numerical geometry within the stated tolerance. Keep its definition, actual checks, measured gaps and repairs. Measure fresh gaps instead of stacking allowances. Apply the class and local limits before admitting a body. For the implemented numerical rules, see [the tolerance code](../crates/catallax-topology/src/numerical/tolerance.rs).
- A bounded answer carries a versioned case rule and an established bound over its full reported scope. A numerical or unproved answer does not become proven or bounded because samples agree. Exact answers, bounded results, numerical results and refusals retain their distinct meanings.
- A refusal explains which operation and condition stopped the request. Never return the unchanged input as success, close a gap silently, or substitute a mesh for a body. Unsupported but possible geometry is unfinished work. A wrong proven answer stops use of that proof method until it is fixed.

## Preview, display and exports

An unproved preview belongs to one immutable recipe revision. The app can replace or confirm it when that revision finishes; late work cannot overwrite a later edit. Keep the preview out of the proof, measurement and export paths. Report a failed build on its own revision and leave an export-specific refusal attached to that request.

An ordinary B-rep part uses its B-rep display path. Field-only features can use a field preview. Render with samples, but use the retained definition for geometry queries and exports. Display meshes share one boundary polyline between incident faces and report a proved error bound where they claim one. Two-kernel sampled agreement is not a whole-boundary proof and does not run on every edit.

Complete rolling-envelope and mitre faces use the shared spline triangulator. It fixes their rail samples before refining face interiors, so neighbouring faces use the same chords. The crossing-saddle strip still uses a separate grid because it collapses at two distinct ends.

Output accuracy is separate from modelling tolerance. STEP preserves exact forms where supported and labels approximate exports. Printing defaults to 0.01 mm; it uses qualified B-rep tessellation for supported solids and certified field output where supported. Keep the requested accuracy and resource limits visible and report a refusal rather than silently coarsening the mesh.

`catallax_app::preview::display_description` reads a retained B-rep into the versioned binary64 input for GPU preview. It returns body-relative carrier definitions, shared edge samples, oriented UV trims and picking names without constructing triangles or running proof. The CPU preview and checked display remain separate calls.

Fitted B-reps use `catallax_app::preview::numerical_display_description`, which reads admitted carriers and retained trim curves and identifies faces and edges by their stable definition keys.

## Work and determinism

Parallel work must not change exact answers, output bytes, topology, bound selection or counted work. Assemble results in a fixed order. Retain deterministic work budgets for interfaces that specify them. The [roadmap](../ROADMAP.md) sets the policy for future app limits and speed measurement; the [verification exchange](VERIFY-INTERFACE.md) describes the separate CLI contract.

Unproved hardware preview values stay outside those proof and export guarantees. IEEE 754 binary64 addition, subtraction, multiplication, division and square root agree under the same operation order, default rounding and subnormal handling, with no intermediate NaN, fused operations or platform transcendental calls. This arithmetic fact does not certify geometry or GPU results.

Kernel crates use safe Rust. Keep third-party types out of kernel interfaces. Adapted code proposes candidates; Catallax checks them and never calls another kernel as a fallback. Keep migrations explicit in the [interface-change guide](INTERFACE-CHANGES.md) and preserve third-party [notices](THIRD-PARTY.md).
