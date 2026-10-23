
Background comparison contract, current coverage and historical whole-check measurements.

## Background comparison prerequisites

`catallax-mesh` now supplies deterministic display tessellation with explicit mesh-to-exact-surface deviation bounds and face picking IDs. `catallax-features` supplies the B-rep recipe adapter for stock definitions and booleans, retaining operation names through trimming and merging. Resolve mesh face IDs through `Model::face_origins` before comparing canonical origin sets; arena IDs are not persistent names. Neither implementation is part of the field crate.

`catallax-compare` supplies the library check. `Controller::begin` invalidates the previous revision and clears its published result. Dispatch `Revision::check` on an application-owned worker, or use `check_built` with a cached B-rep build and the field recipe. Only a complete report for the still-current revision can publish. `cancel` invalidates its own revision without touching a newer edit. No library thread pool, clock deadline or app-side edit tracking is added.

The forward pass evaluates the certified field query at every selected mesh vertex, preserving exact binary64 coordinates as rationals. If a residual supplies no finite upper, prepared field-native trimmed-face queries seek exact admitted witnesses. For mesh deviation `e` and field distance interval `[l,u]`, the comparison uses the lower bound `max(0,l-e)` and upper bound `u+e`. A lower bound above the display accuracy proves a gap. An exact topological endpoint may supply another candidate, but only after the independent field proves that point is on its result boundary; the measured point-to-witness distance then supplies a certified upper bound. An unbounded upper remains unresolved when no such witness settles it. A supporting-surface hit without a surviving trim witness never proves proximity, and no raw Boolean residual decides it.

Incident B-rep origins guide field witness order, with curved supports first.
They never establish trim membership or exclude other field faces. A failed
preferred witness falls back to the complete field registry.

The reverse pass generates rational stock-boundary and feature-boundary candidates independently of the B-rep mesh. `Plan::visit_feature_candidates` uses the field's stored profile and sweep definitions, including signed and oblique extrusion, rational-frame full/partial revolution and through-all walls. Through-all sampling uses a conservative field-result bounding box propagated through the expression tree, not artificial tool caps or a B-rep bound. Certified field face queries discard removed candidates and retain unresolved trimming. Each surviving sample first tries an exact rational support witness and checks it against the matching trimmed B-rep face. If that projection misses the trim, `catallax_measure::distance::PreparedFace` refines certified analytic chart boxes and admits exact rational or algebraic witnesses through the boundary-loop predicate. Lower certificates cover every retained cell; upper certificates use only admitted face points. These queries change neither kernel's generated geometry.

Origin comparisons use canonical sets from `Model::face_origins`. A maximal merged face can contain different creating origins on adjacent patches, so the check collects the locally witnessed names over that face rather than requiring every name at every point. Field normal enclosures and exact B-rep gradients check outward direction at the reverse witnesses. Sharp mesh vertices do not acquire fabricated unique normals. Missing origin or direction coverage stays unresolved.

`Report::sampled_distance_upper` is the maximum checked distance upper, including forward mesh error. `None` means at least one checked distance remains unbounded. The verification CLI publishes it as `comparison.sampled_distance_upper_mm`. It is a sampled bound, not a whole-part Hausdorff certificate. See [trimmed minimum distances](FIELD.md#trimmed-minimum-distances) for the shared interval engine and the independent trim adapters.

Defaults are 0.01 mm display accuracy, 0.005 mm mesh chord tolerance, four rational intervals per sampling parameter interval, and minimum normal cosine 999/1000. These are sampled display checks, not a global Hausdorff certificate or the finer 1 nm pre-export check. Raising `divisions` increases independent reverse sampling. A closed rational `Region` restricts both passes; forward selection pads the box by the local mesh error. A scope with no two-way samples stays unresolved. A merged-face origin whose extent cannot be established within the region stays unresolved rather than causing a remote naming change to fail the local check. The library still rebuilds and meshes the complete B-rep, so regional rechecks reduce query work, not build cost. Cancellation is cooperative between kernel calls and individual samples; existing build and meshing calls are not interrupted internally.

| Recipe operations | Comparison coverage |
|---|---|
| Block, cylinder, sphere, cone frustum, ring torus | Both builds, display vertices, independent exact field samples, trimmed-face projections, names and directions |
| Unite, subtract, intersect | Same two-way path wherever both kernels build and mesh the recipe; unresolved Boolean query bounds remain unresolved |
| Sketch, extrude, extrude-through, revolve | Same two-way path for the shared rational support subset, including profile holes and segment aliases; specific unfinished frame/surface refusals stay unresolved |
| Linear/circular body patterns and mirrors | Same two-way path, including nested copy origins and transformed stock/feature candidates; rational-placement curved copies use independently admitted trimmed-face witnesses, while algebraic-only sampling and numeric refusals remain unresolved |
| Path, hole, fillet, chamfer, shell, offset, thicken, draft, sweep, loft, STEP import | Not comparable: neither recipe adapter implements the operation |
| Lattice, gyroid, smooth blend, mesh import | Not comparable: the field implementation is missing; a B-rep implementation is not planned |

Every unsupported step is returned with its operation-specific reason, not only the first step. Builds, meshes and query proofs that fail stay unresolved. A proven discrepancy takes precedence over unresolved observations, which remain in the report. Exact expression identities and opposite coincident primitive halfspaces can certify agreement of two empty results without boundary samples. General emptiness is not inferred from missing samples.

Named field queries are fallible even when their handles are valid. A numeric refusal in reverse sampling records its code as an unresolved observation and stops that sampling phase; it is never an internal panic or evidence that a face disappeared. The verification CLI discards a request whose shared work allowance is exhausted and returns `numbers.operation-budget`, rather than publishing partial comparison results.

### Coplanar booleans and repeated operands

The joined-blocks recipe below now agrees with no unresolved observation. The field proves the removed meeting plane interior rather than treating its zero residual as boundary.

The coplanar reproducer uses two size-2 blocks at translations `[0,0,0]` and `[1,0,0]`, with identity rotations. Union, subtraction and intersection all agree without unresolved observations. For subtraction, samples such as `[3/2,0,1/2]` are outside the regularised result and are not selectable source faces. Flush blocks translated by `[2,0,0]` also agree for all three operations, including their empty intersection.

A Boolean may name the same body step twice. Union and intersection retain its material and normalise nonmaximal topology; subtraction produces an empty body. The block `a` of size `[2,2,2]` at identity placement, followed by subtracting `a` from `a`, builds and compares as empty without a capability exclusion.

Remaining meshing limitation:

- The default one-million-vertex meshing allowance is insufficient for a torus with major radius 30 and minor radius 10 at 0.005 mm chord tolerance. The failure is `mesh.precision-or-budget`, not impossible geometry or disagreement. Reproducer: stock torus `stock`, identity rotation and translation `[7/3,-5/7,2]`. Periodic temporary chart copies count toward this allowance. The agreeing scale fixtures request four million vertices without changing the 0.01 mm comparison accuracy or 0.005 mm mesh tolerance; production defaults are unchanged.

The application still owns edit-region calculation, worker scheduling, restart policy, attribution to the edited feature and the kernel bug list. SolidDraft integration, development checks, pre-export checks and command-line exchange changes are outside this library.

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
These historical measurements predate the coplanar Boolean fixes above. In
particular, the joined-blocks outcome in this table is not its current outcome.

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

No unperturbed representative part produced a proven discrepancy. At measurement
time, joined blocks exposed indeterminate field queries on the removed meeting plane.
At that time, curved Boolean mesh vertices generally had unbounded field distance
uppers; those observations did not become agreement. Current limitation regressions
instead require finite trimmed bounds, while changed dimensions, names and
orientation still require proven discrepancies.

No general interactive budget or candidate-versus-baseline slowdown is established
by these shared-machine samples. Attribute a performance regression only after
three matched alternating baseline/candidate pairs on one pinned core, recording
load and geometric work units. Finite trimmed-face upper bounds now provide the
distance proof missing from the historical curved examples. Application integration
still needs its own edit-loop and interactive latency evidence.
