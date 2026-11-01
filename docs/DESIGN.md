# Catallax design rules

[PRD.md](PRD.md) sets the product direction. This document governs the B-rep half of Catallax. [ARCHITECTURE.md](ARCHITECTURE.md) describes its implementation, and [pieces/](pieces/) records the existing crate contracts. [PLAN.md](PLAN.md) orders the work still needed.

## Fixed points

1. One recipe, the feature tree, is the master. An exact B-rep kernel and a signed distance field kernel build each part from that recipe. Every operation has one written mathematical definition shared by both kernels.
2. The B-rep stores shapes, edges and vertices by their definitions. Coordinates and approximating curves are caches, never geometric truth. A proven answer is never wrong.
3. The kernels agree exactly where mathematics allows. Where necessity forces approximation, the proven geometric discrepancy is at most 1 nm, or 0.000001 mm. Section 3.1 of the PRD defines the only two permitted cases.

The repository implements the B-rep starting code and an independent analytic field evaluator for recipe stock solids and booleans. [RECIPE.md](RECIPE.md) defines the shared language, `catallax-recipe` holds its data model, and [FIELD.md](FIELD.md) states the field query guarantees and unresolved cases. B-rep recipe evaluation, agreement checks and operation-derived B-rep names are still planned. The existing B-rep numerical and exchange contracts are unchanged.

## Rules for the inside of the B-rep kernel

### R1. Every decision is certified or marked

A decision changes what the kernel builds: a sign, side, intersection or ordering question. Every decision goes through piece 1 and has one of four outcomes:

- **Proven.** The strict answer is certified.
- **Proven coincident.** Exact zero has a proof.
- **Not exact, with a proven bound.** An unfinished proof may take the coincident branch only after proving the resulting geometric gap is smaller than 1 nm. The mark travels with the result. A raw scalar or angle bound is not a distance certificate; it must be converted and propagated into a bound on the affected geometry.
- **Named refusal.** Stop without publishing a partial body. Only impossible geometry is an acceptable final refusal. Any refusal of possible geometry is unfinished implementation work.

Examples of impossible input include a non-positive radius, a self-intersecting torus and a connected shell touching itself along a curve. Unsupported pairs, ambiguous imprints, exhausted budgets, precision limits, timeouts and unproven measurement accuracy are bugs or missing coverage, not acceptable exclusions from the product.

No code outside piece 1 compares floating-point numbers to decide topology. Floating-point values may occur in error-bounded filters and caches for display, export and search. No guessed tolerance decides a topological question.

The inherited APIs can still return bounds larger than the new 1 nm requirement, including non-distance marks. Enforcing and composing that geometric limit across all operations is open work. A marked result today is not a claim that the new requirement is already met.

### R2. A proven answer is never wrong

A wrong proven or proven-coincident result stops work on the responsible proof method until the cause is repaired and checked on fresh mathematical cases. Refusals and marked answers are reported separately from wrong proofs. They are not evidence of completed feature coverage.

### R3. Work budgets are deterministic

Use the fast certified filter first, then increased precision and exact methods. Count budgets in deterministic work units, not elapsed time. Exhaustion permits a bounded outcome only under R1's geometric limit; otherwise it exposes unfinished work. Every accepted mark is visible, counted and carried on the result.

**Open conflict, no behaviour change in this publication.** The CLI still sets a 50-second wall-clock deadline and reports `numbers.deadline`. Machine load can make the same question answer on one run and refuse on another. This conflicts with R10's same-answer rule. The existing exchange contract documents the exception; resolving it needs a separate design and versioned contract decision.

### R4. Store by definition, coordinates as a cache

Planes, cylinders, cones, spheres, tori, lines, circles and ellipses have exact rational parameters. A derived intersection curve records its supporting surfaces and branch witness. A vertex records its defining meeting and which solution it is. A future fillet records the rolling ball and its support faces; an offset records the source and distance.

Coordinates, parameters, boxes and fitted splines may be discarded and recomputed at any precision. Decisions read the definition or a proven enclosure, never an unverified cached coordinate.

### R5. Record construction facts

Record tangency, coplanarity, coaxiality, concentricity, parallelism and incidence only when guaranteed by construction or exact proof. Consult facts before numeric refinement. Never infer a fact from a near-zero value or a marked outcome.

### R6. Curved topology remains general

Topology never assumes planar faces, straight edges, three-plane vertices or a distinguished outer loop. A face can have no loops, several equal-standing loops, or a point loop at a cone apex. An edge can be a closed ring without vertices. Leave room for later field features without pretending they are already implemented body kinds.

### R7. One boolean engine

Intersect, imprint, classify, keep and sew are the five shared stages. Later topology-changing features reuse them rather than growing separate boolean implementations.

### R8. No angles inside decisions

Use exact cyclic-order predicates for circles, ellipses and closed branches. Transcendental parameters are caches for display, export and measurement, not decision inputs.

### R9. Scale and placement do not change answers

Exact rational scaling and exact rigid moves preserve topology and transform the answer accordingly. Test scaled and moved copies. The 1 nm allowance is a ceiling on necessary approximation, not permission to snap geometry at an absolute size threshold.

### R10. Determinism

The same recipe and question must produce the same outcome, topology and proven bounds on every machine and at every thread count. Do not let hash iteration, scheduling or platform `libm` choose results. Use the kernel's own proven transcendental bounds. The current deadline conflict remains open under R3.

### R11. Refuse loudly, never silently

A current failure reports its piece, stage and stable code. Never silently return the input, fill a gap, move geometry or substitute a mesh. Impossible geometry is a legitimate refusal; every other refusal is unfinished work. Bad external input is a structured failure, never a panic. Invalid handles and violated internal preconditions are programming errors.

### R12. Speed per click matters

Catallax serves an interactive CAD tool. Measure rebuild and query latency on representative edits, including refusal latency, and reduce costs that block a click. Correct proofs take precedence over speed. Timing measurements are not correctness thresholds, and no current test fails just for being slow. Background rebuilds and safe parallel execution are planned; neither may change answers.

### R13. Safe Rust

Kernel crates forbid `unsafe`. Growing public enums are non-exhaustive. Dependencies point down the crate graph; no external dependency type crosses a public boundary without review.

### R14. Name entities by their creating operation

A future recipe reference names, for example, the edge between the side of extrude 1 and the floor of cut 2. It must not depend on result traversal order or arena index. Existing store-local handles remain implementation identities, not persistent feature references. [RECIPE.md](RECIPE.md) section 3 defines the naming and rebuild rules; the B-rep does not report these names yet.

## Conversion at the border

Conversion is built only when an operation or output needs it, not after every edit. Imported bodies become fixed recipe inputs. The [STEP importer](pieces/08-step-import-repair.md) reads a supported analytic subset and measures and records straight planar seam repairs. Self-exported blocks and cylinders pass through the application recipe entry with unchanged carriers and boundaries. General periodic faces, curved repairs and spline import remain unfinished. Any movement must be measured against the import's tolerance and reported.

For export, keep analytic forms exact. A surface with no exact STEP form stays exact by definition in both kernels; only its exported representation is approximated. Prove and report fitting and rounding deviation, subject to the PRD's 1 nm necessity limit. Existing configurable export tolerances and marks are an inherited contract, not proof that every exported result already meets this new rule. Accuracy levels and large-part stand-ins remain planned.

## Target and evidence

The target is reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels. There is no claim to beat another kernel. The former 1-in-10,000 programme and grid-cell acceptance scheme are retired as product goals. Existing generated and degenerate tests remain useful regression evidence, not statistical product guarantees.

Keep exact mathematical cases, independent oracles, scale and placement checks, repeatability checks and measured timings. An independent tester built without seeing the kernel's code checks the public exchange contract. No result here claims complete feature coverage or parity with another kernel.

## Changing contracts

Propose public API and format changes through review. Update the affected piece descriptions and add a new [interface-change fragment](interface-changes/README.md) with migrations for every caller. Preserve the legacy exchange names until a versioned contract change replaces them. New PRD requirements must identify inherited implementation gaps rather than relabel the existing code as compliant.
