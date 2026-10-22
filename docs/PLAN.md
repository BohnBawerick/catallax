# Catallax implementation direction

[PRD.md](PRD.md) defines the product. This plan separates the inherited B-rep implementation from the work still required. The target is reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels. The former 1-in-10,000 programme and grid-cell acceptance goals are retired; keep their useful mathematical regressions.

## Keep the B-rep starting code

Carry forward exact arithmetic and certified decisions with deterministic work budgets, definition-based analytic geometry, append-only curved topology, the five-stage boolean engine, bounded mass properties and STEP AP242 export. Preserve the public exchange contract unchanged, including its legacy names. `catallax-features` adds named B-rep recipe stock solids and rationally represented sketch features. The independent `catallax-field` evaluator covers stock solids and booleans; see [FIELD.md](FIELD.md) for its guarantees.

The existing tests include ignored cases. Their source records the reasons; a green default suite does not mean those cases were repaired. Existing work still in the source repository must arrive as a separately reviewed export delta, never by silently replacing this tree.

## Next work, in order

1. Continue independent testing of the B-rep half and repair refusals of possible geometry. General cone and torus pairs, ambiguous imprints and classification, precision limits and missed measurement accuracy remain unfinished work. Only impossible geometry is an acceptable final refusal.
2. [RECIPE.md](RECIPE.md) defines each operation and its origin-based naming. Named B-rep recipe builds now preserve origins through splits and maximal merges. Settle the remaining section 7 items and widen exact geometry beyond rational analytic frames. Arena indices are not persistent names.
3. Define and enforce a composable geometric error policy. A bounded decision may stand in for exact coincidence only after proving the gap smaller than 1 nm. Existing scalar, angle and distance marks do not yet establish that invariant over a whole result.
4. The field evaluator covers stock solids, sketch profiles, signed extrusions, two-sided through-all cuts, rational-frame revolutions and booleans. `catallax-compare` performs sampled two-way display comparison over that common subset. Extend certified classification of degenerate regularised zeros and finite closest-point bounds on surviving trimmed faces before claiming complete distance/proximity queries or generally resolved curved-sweep and Boolean comparisons.
5. Finish the remaining sketch, extrusion and revolution geometry bugs listed in RECIPE.md section 7.14 in both evaluators, preserving common definitions and comparable origin sets. Continue through holes, fillets and chamfers, shell and offset, draft, sweep and loft, pattern and mirror. Every operation needs its shared definition before two implementations can agree.
6. Integrate the shared-edge display meshes from `catallax-mesh`, then add STEP import, on-demand conversions, output hand-off checks, export accuracy levels and large-part stand-ins. Seam curves require a versioned exchange-contract change. Measure export interoperability rather than assuming a syntactically valid file works in every CAD tool.
7. Integrate the kernel into SolidDraft once the required operations are covered. Add field-only features, printing, FEA queries and interference checks with the hand-off guarantees in the PRD. FEA solver design remains separate.

## Feature build order after the shared foundation

RECIPE.md sections 7.1 to 7.13 now record the settled definitions and product choices, with future variable-radius and nonplanar parting-line extensions kept distinct. Its version 4 specification is docs-only; the crate still reads version 3. Land the schema and migration review separately before emitting version 4 files. Do not relabel arc-length lofts automatically.

1. Finish algebraic sketches, extrusion and revolution, then named holes and the additional extrusion end conditions. Root-defined sketch values, frames, rational-degree angles, elliptic cylinders and analytic trims consume the shared `Real`/`RealPoly` foundation. First-hit caps need certified trimmed-face intersections and finite bounds. Pattern and mirror can follow as soon as their inputs and transformed carriers work; saved feature-pattern expansion needs no separate geometric kernel.
2. Add complete trimmed-face distance queries and regularised offset construction. Build distance shells, outward-shell compositions and thicken from those queries, first on analytic reductions, then on certified pipe and envelope carriers. Keep straight thicken rims and shared origin sets. A sign-only field residual cannot supply these metrics.
3. Build the defined analytic fillet/chamfer runs and planar/cylindrical draft cases. Then implement the agreed roll-over/keep-neighbours modes, corner continuation, consumed-face history and neutral-plane draft retention with kept-region overrides. Do not claim these paths from support equations alone.
4. Land generic derived carriers with their charts, predicates, intersections, integration and export together. Then add rational lofts, fixed spline sweeps and general algebraic offset/blend/draft supports. Analytically reducible line/arc sweeps and their stationary round turns can land earlier. Helices and general rotation-minimising transport need an additional certified function/ODE representation beyond algebraic `Real`.
5. Implement STEP replacement edits and persistent-reference editing alongside import support. Add reference-only open scans and reconstruction deviation reports alongside closed mesh import; never infer a solid from the scan. No edit may infer identity from a distance tolerance alone.

At each stage compare point sets, surviving origins and refusal classifications in both kernels. Require whole-result error transfer proofs under ARCHITECTURE.md before accepting bounded geometry; display agreement is separate evidence. Reassess after shared analytic operations and their comparison proofs land, before committing to general derived carriers.

## Background comparison prerequisites

`catallax-mesh` now supplies deterministic display tessellation with explicit mesh-to-exact-surface deviation bounds and face picking IDs. `catallax-features` supplies the B-rep recipe adapter for stock definitions and booleans, retaining operation names through trimming and merging. Resolve mesh face IDs through `Model::face_origins` before comparing canonical origin sets; arena IDs are not persistent names. Neither implementation is part of the field crate.

`catallax-compare` supplies the library check. `Controller::begin` invalidates the previous revision and clears its published result. Dispatch `Revision::check` on an application-owned worker, or use `check_built` with a cached B-rep build and the field recipe. Only a complete report for the still-current revision can publish. `cancel` invalidates its own revision without touching a newer edit. No library thread pool, clock deadline or app-side edit tracking is added.

The forward pass evaluates the certified field query at every selected mesh vertex, preserving exact binary64 coordinates as rationals. For mesh deviation `e` and field distance interval `[l,u]`, the comparison uses the lower bound `max(0,l-e)` and upper bound `u+e`. A lower bound above the display accuracy proves a gap; an unbounded upper bound remains unresolved even when an auxiliary surface projection exists. No raw Boolean residual decides proximity.

The reverse pass generates rational stock-boundary and feature-boundary candidates independently of the B-rep mesh. `Plan::visit_feature_candidates` uses the field's stored profile and sweep definitions, including signed and oblique extrusion, rational-frame full/partial revolution and through-all walls. Through-all sampling uses a conservative field-result bounding box propagated through the expression tree, not artificial tool caps or a B-rep bound. Certified field face queries discard removed candidates and retain unresolved trimming. Each surviving sample projects to an exact rational supporting-surface witness, then checks that witness against the matching trimmed B-rep face. Analytic support distances supply lower bounds; the measured distance to a certified on-face witness supplies the upper bound. A chart projection need not be the closest point. The added queries do not change either kernel's generated geometry.

Origin comparisons use canonical sets from `Model::face_origins`. A maximal merged face can contain different creating origins on adjacent patches, so the check collects the locally witnessed names over that face rather than requiring every name at every point. Field normal enclosures and exact B-rep gradients check outward direction at the reverse witnesses. Sharp mesh vertices do not acquire fabricated unique normals. Missing origin or direction coverage stays unresolved.

Defaults are 0.01 mm display accuracy, 0.005 mm mesh chord tolerance, four rational intervals per sampling parameter interval, and minimum normal cosine 999/1000. These are sampled display checks, not a global Hausdorff certificate or the finer 1 nm pre-export check. Raising `divisions` increases independent reverse sampling. A closed rational `Region` restricts both passes; forward selection pads the box by the local mesh error. A scope with no two-way samples stays unresolved. A merged-face origin whose extent cannot be established within the region stays unresolved rather than causing a remote naming change to fail the local check. The library still rebuilds and meshes the complete B-rep, so regional rechecks reduce query work, not build cost. Cancellation is cooperative between kernel calls and individual samples; existing build and meshing calls are not interrupted internally.

| Recipe operations | Comparison coverage |
|---|---|
| Block, cylinder, sphere, cone frustum, ring torus | Both builds, display vertices, independent exact field samples, trimmed-face projections, names and directions |
| Unite, subtract, intersect | Same two-way path wherever both kernels build and mesh the recipe; unresolved Boolean query bounds remain unresolved |
| Sketch, extrude, extrude-through, revolve | Same two-way path for the shared rational support subset, including profile holes and segment aliases; specific unfinished frame/surface refusals stay unresolved |
| Path, hole, fillet, chamfer, shell, offset, thicken, draft, sweep, loft, pattern, mirror, STEP import | Not comparable: neither recipe adapter implements the operation |
| Lattice, gyroid, smooth blend, mesh import | Not comparable: the field implementation is missing; a B-rep implementation is not planned |

Every unsupported step is returned with its operation-specific reason, not only the first step. Builds, meshes and query proofs that fail stay unresolved. A proven discrepancy takes precedence over unresolved observations, which remain in the report. Exact expression identities can certify agreement of two empty results without boundary samples. General emptiness is not inferred from missing samples.

Existing kernel limitations exposed during comparison development remain unrepaired here:

- Coplanar subtraction can return indeterminate field trim queries on removed common supports. Reproducer: two size-2 blocks at translations `[0,0,0]` and `[1,0,0]`, with identity rotations, then subtract the second from the first. Samples such as `[3/2,0,1/2]` remain unresolved, not disagreements or agreement certificates.
- A Boolean naming the same body step as both operands reaches the B-rep engine's distinct-body assertion. The comparison capability report marks that case not comparable before entering the assertion. Reproducer: block `a` of size `[2,2,2]` at identity placement, then subtract `a` from `a`. Two separately named identical blocks do build and compare as an empty result.
- The default one-million-vertex meshing allowance is insufficient for a torus with major radius 30 and minor radius 10 at 0.005 mm chord tolerance. The failure is `mesh.precision-or-budget`, not impossible geometry or disagreement. Reproducer: stock torus `stock`, identity rotation and translation `[7/3,-5/7,2]`. Periodic temporary chart copies count toward this allowance. The agreeing scale fixtures request four million vertices without changing the 0.01 mm comparison accuracy or 0.005 mm mesh tolerance; production defaults are unchanged.

The application still owns edit-region calculation, worker scheduling, restart policy, attribution to the edited feature and the kernel bug list. SolidDraft integration, development checks, pre-export checks and command-line exchange changes are outside this library.


## Deterministic verification policy

The policy conflict is closed for the opt-in `catallax-verify` exchange. Whole-request work budgets choose outcomes; wall-clock measurements report performance only. External cancellation leaves an incomplete run. The legacy CLI retains its 50-second deadline unchanged until callers migrate. [VERIFY-INTERFACE.md](VERIFY-INTERFACE.md) defines the versions, adapters and replay driver; [ARCHITECTURE.md](ARCHITECTURE.md#the-time-budget) records the allowance and interactive targets.

## Measure useful latency

Record rebuild and query time per representative edit, including refusal latency and agreement-check cost. Speed per click matters. Optimize measured costs without weakening proofs. Safe background threads and cancellation must preserve deterministic answers and must never expose a half-built result to the application.

### Whole-check measurements, 2026-10-22

Measured on an x86_64 Linux workstation with an Intel Core i7-9750H, release
build, process priority `nice 10`, default comparison options. Each row is one
complete run including independent field compilation, named B-rep construction,
display meshing, both query passes and revision-gated publication. Compilation
time is excluded. Other validation work was running concurrently; these are
observed shared-machine costs, not idle-machine throughput estimates.

```sh
nice -n 10 cargo run --release -p catallax-compare --example measure -- 1
```

The example fixes the complete recipes, including identity placements. Dimensions
are millimetres. The joined blocks share a complete end face. The bored block
subtracts a radius-1.5 cylinder of height 7 at `[5,5,-1]`. The intersecting spheres
have centres `[0,0,0]` and `[2,0,0]`.
The stock measurements preceded the field-feature integration rebase. The three
profile rows were measured afterwards using the corresponding example case
filter. All rows use the same display accuracy and meshing controls.

| Example | Dimensions | Whole check, ms | Forward vertices | Reverse samples | Outcome | Unresolved observations |
|---|---|---:|---:|---:|---|---:|
| `block` | 20 x 10 x 5 | 20.145 | 8 | 96 | Agree | 0 |
| `cylinder` | radius 2, height 5 | 1393.344 | 496 | 192 | Agree | 0 |
| `sphere` | radius 2 | 19674.712 | 8066 | 128 | Agree | 0 |
| `frustum` | radii 2 and 1, height 4 | 11564.849 | 4376 | 192 | Agree | 0 |
| `torus` | major 3, minor 1 | 52706.842 | 16384 | 256 | Agree | 0 |
| `joined_blocks` | 20 x 10 x 5 united with 5 x 10 x 5 | 275.867 | 8 | 160 | Unresolved | 32 |
| `bored_block` | 20 x 10 x 5 with a through bore | 2163.574 | 400 | 128 | Unresolved | 178 |
| `intersecting_spheres` | two radius-2 spheres | 51719.860 | 3294 | 44 | Unresolved | 3293 |
| `profile_extrude` | 4 x 3 rectangular profile, signed height 2 | 38.059 | 8 | 96 | Agree | 0 |
| `profile_through_slot` | 20 x 10 x 5 block, 4 x 2 rectangular through slot | 211.469 | 16 | 158 | Agree | 0 |
| `profile_quarter_revolve` | radii 1 to 2, height 3, 90-degree revolution | 1041.295 | 208 | 96 | Unresolved | 106 |

No unperturbed representative part produced a proven discrepancy. The joined
blocks expose indeterminate field queries on the removed meeting plane.
Curved Boolean mesh vertices generally have unbounded field distance uppers;
those observations do not become agreement. The tests also deliberately change
dimensions, locations, names and orientation to require disagreements.

No general interactive budget or candidate-versus-baseline slowdown is established
by these shared-machine samples. Attribute a performance regression only after
three matched baseline/candidate rounds on separate pinned physical cores,
recording load and geometric work units. Field Boolean queries also need finite
upper distance bounds before the background check can routinely settle those
edits. That is a useful reassessment point before app integration.

## Evidence for completion

Use mathematical cases, independent numeric oracles, scale and placement variations, repeatability checks and measured timings. Keep refusals, marked answers and wrong proofs distinct. An independent tester built without seeing the kernel's code uses the documented exchange contract. No unpublished tester design or acceptance system forms part of this plan.
