# M19 bounded adaptive print request and sheet fields

Baseline: `39aa85a8513ca1b0283a5f28c188b74df6d08412` (current main fetched
11 November 2026). **This is a partial M19 increment, not milestone completion.**

## Sheet definition

An explicit `kind: "sheet"` now compiles to the regularized difference of two
existing network fields, `N(c) −* N(−c)`, clipped to the requested region. This
is exactly the band `|g| <= c`, not a thickness in millimetres. It reuses the
existing exact point/ray decisions, polynomial interval enclosures and Boolean
semantics; there is no new sampled membership rule or trigonometric evaluator.

Constant levels must be positive. Linear grades must have positive endpoints
and retain the existing network restrictions: one-axis anchors, endpoints at
most 1/4 and slope at most 1/10 per millimetre. Negating both endpoints commutes
with clamping. Nonpositive sheet levels refuse with `gyroid.impossible`;
unsupported frames/grades keep typed coverage refusals. A constant band covering
the whole network range reduces to its region.

Terminal stage construction accepts these sheet definitions. **Sheet printing
remains refused** with `field.print-detail-topology` in legacy print routes:
the existing single-network certificate does not establish both local boundary
sheets. A private provenance flag marks only live leaves of explicit sheet
operations; pre-existing network Boolean recipes retain their previous route.
The sheet implementation does not add a mass/stand-in solver or qualify those
outputs. Existing consumers may refuse their unsupported two-network graph.

## Explicit adaptive request

`printing::adaptive::mesh(plan, options, dec)` and `mesh_with_workers` return
an immutable `Output` or `Box<Failure>`. Existing `printing::Options` retain
their meaning and default 0.01 mm accuracy. Existing printing entry points
and network/lattice outputs are not redirected to this new route.

The candidate uses the landed conforming rectangular grid and Kuhn tetrahedra.
Refinement splits the longest implicated cell axis; a split propagates as an
entire coordinate plane. Thus neighboring cells share identical face vertices
and diagonals even where widths differ; there are no hanging nodes to weld by
a tolerance. Plane propagation can refine much more than a balanced octree;
this is a bounded first cut, not an optimal adaptive mesh.

The new request shifts its initial candidate by 1/32 of the coarse pitch to
avoid common exact clip/grid contacts. This is only candidate selection, not
a tolerance change. Certification still uses the original requested accuracy.
Unknown cells with homogeneous corner signs are retained for refinement; they
are not inferred empty and cannot disappear because a neighboring facet is
close. The existing three-refinement limit, coordinate constraints and cell/
vertex caps remain. Failure publishes no mesh.

The inherited two independent interval proofs cover every mesh patch to the
true boundary and the complete true boundary to the mesh. Only after both
succeed does each emitted cell retain exact bounds, a contiguous facet range
and the conservative certified error. This bound is presently the same global
bound for each cell; it is not a tighter estimate derived from cell size.
For true boundary points in a cell it bounds distance to the **complete mesh**,
not necessarily to that cell's own facets. Strictly classified interior/
exterior cells emit no facets. The report is bound to immutable geometry.

**General material-topology certification remains a blocker.** Hausdorff
coverage and a watertight candidate alone can lose a tiny void inside a mixed
cell. The new API therefore admits only a single convex block (or proven empty),
requiring one genus-zero mesh component. Other otherwise certified candidates
refuse at `MaterialTopology`; unresolved homogeneous cells can refuse earlier.
The implementation does not assert that this guard proves general lattice or
sheet printability. P6, thin-feature and sheet acceptance remain open, rather
than being silently weakened to a distance-only test.

## Refusal and concurrency contract

`Failure` retains the original `Refusal`, including attached undecided evidence,
the unchanged request, the active stage and this call's Decider arithmetic work.
Stages distinguish preparation, candidate topology, mesh-to-boundary coverage,
boundary-to-mesh coverage, refinement, material topology and representation.
Decider work is not total request work: `request_work` also retains the charged
units including cancellation checkpoints when a caller's operation allowance is
active, and is explicitly unavailable otherwise.

A positive-accuracy failure offers a **separate UNVALIDATED coarser request** at
twice the tolerance with unchanged caps. It is not executed automatically and
is explicitly **not estimated to fit**. Cells, facets, peak memory and calibrated
time estimates remain unavailable; no throughput calibration was performed in
this shared environment. A coarser request cannot solve missing topology proofs.
There is no claimed higher-cap price or fabricated timing prediction.

Default worker count is `available_parallelism` (four on this executor).
Independent final facet representation checks actually run on processor workers.
Fixed batches of at most 1024 facets are prepaid on the caller's Decider; worker
code uses only fixed-cost output-coordinate arithmetic, not numeric Deciders or
thread-local request allowances. Every worker is joined and cancellation is
checked before publication. Geometry, certificate order, failure precedence and
charged work are independent of worker count. Thread-start/panic failures are
typed infrastructure failures, not geometric conclusions.

**Meshing and interval certification remain serial.** General parallel numeric
proof work cannot safely inherit the current thread-local allowance implicitly.
No fresh per-worker proof budget was introduced. The all-core final validation
stage is not evidence of all-core meshing speed-up; that assignment exit remains
open. Tests at finite thread counts are not universal scheduling or cross-platform
proofs, and wall-clock deadlines may naturally expire at different stages.

## Validation and remaining exits

The public `adaptive_print` Rust controls cover exact positive/negative band
membership and both sheet boundaries, axis grades, impossible levels, unequal
cell widths, complete facet ownership and local bounds, an independent M22
closed-mesh/material check, worker equality of full geometry and certificates,
exhausted work allowances and expired deadlines, stage/choice labels, tiny-cavity
refusal and legacy sheet export refusal. Refusal equality includes attached
undecided evidence, not just Debug strings.

Public P6 and V7c network requests also exercise the new refusal report at the
unchanged 0.01 mm accuracy with a deliberately small 256-cell cap. Those are
capacity controls, **not full P6/gyroid qualification**. Existing network and
legacy printing focused regressions are separate preservation checks.

Observed on 11 November 2026 with Rust 1.99.1, x86_64 Linux, release profile:

| Control | Result | Charged request work |
|---|---|---:|
| 4 × 1 × 1 mm block, 0.125 mm request, 1/2/4/default workers | Same 4,368 facets, 13,710 visited cells and complete cell evidence; independent M22 validation passed | 63,799 each |
| Same block, one unit below successful work | Same representation-stage refusal at every tested worker count | Allowance exhausted |
| Sheet control, period 4, level 0.25 | Refinement refusal: 436 unresolved homogeneous cells | 693,843 |
| Radius 0.03125 cavity control | Refinement refusal: 5 unresolved homogeneous cells | 97,015 |
| P6 and V7c, 0.01 mm, 256-cell cap, 1/2/4/default workers | Same preparation-stage `field.print-budget`; accuracy/caps unchanged | 1 each |

Six new tests, 15 existing printing tests and three existing gyroid unit tests
passed. Package all-target Clippy, workspace formatting check and whitespace
check passed. The new tests also compare complete failure values, including
original request and undecided evidence, under 1/100/1000-unit allowances.
No elapsed-time number is used as calibration or a gate. Raw command logs are
retained outside source in the saved execution environment.

Reproduction:

```sh
cargo test --release -p catallax-field --test adaptive_print -- --nocapture
cargo test --release -p catallax-field --test printing
cargo test --release -p catallax-field --lib gyroid::tests
cargo clippy -p catallax-field --all-targets -- -D warnings
cargo fmt --all --check
git diff --check
```

Full P6 at its recorded cap/accuracy, printable sheet output, general mixed-cell
material equivalence, parallel interval/refinement work, calibrated trade-off
estimates, and a fit-certified coarser option remain unestablished. No 2-second
timing gate applies to these results. No C++/OCCT, private/held-out inputs,
saved-answer inventories, allowance increases or full-workspace/CI reruns were
used. Hosted final-head integration and independent review remain separate.
