# F7 no-mesh solver queries

This is an opt-in Rust interface in `catallax_field::solver`, built on the
merged [FEA design](FEA-DESIGN.md) and F6 classifier. It supplies geometry and
finite-cell integration data, not certified stress or a production elasticity
solver. M16 and M20 are not complete. No proof/export route consumes solver
estimates, and no existing tolerance or operation allowance changes.

`Definition::compile(recipe, materials, decider)` owns the native field plan.
Its revision hashes the canonical recipe and length-delimited material identity
and coordinates. `prepare` produces shared immutable `Queries`; callers retain
private `Workspace`, `BoxStorage` and `Decider` objects. No B-rep is needed by the
adapter. Imported file maps are not accepted by this constructor. Recipe and
numeric coverage failures remain typed refusals.

All coordinates/distances are mm, volume/moments are mm³ in normalized local
coordinates, area is mm², forces are N and uniform tractions are N/mm². Materials
are explicit identities, not silently supplied densities or elastic constants.

## The six services

| Service | Contract and coverage | Uncertainty/refusal |
| --- | --- | --- |
| Q1 `membership` | Delegates exact rational points to `Plan::proven_checked`; keeps the distinction between residual and true distance. | Boundary and Indeterminate remain separate from Inside/Outside. Numeric refusals retain their attached evidence. |
| Q2 `whole_cell` | F6 proof over the entire closed rational axis box. | Unknown is unresolved, not a proved cut; zero/negative cell extents refuse. |
| Q3 `integration`, `Cell::moments` | Exact clipping for a single native axis box with signed axis permutation/translation. Otherwise use F6 on the cell, then eight octants of an unknown cell. Return attributable material pieces and fraction bounds. All 27 tensor moments are available per region. | Unresolved pieces contribute full upper measure and zero lower measure. No requested integration width is implied. Connectivity is certified only for a nonempty native-box intersection. General cut-cell stiffness must refuse unresolved connectivity/integration quality. |
| Q4 `box_distance`, `boundary_distance`, `distances().point_to_boundary` | Native box: exact squared true-boundary distance and all closest displacement vectors, with an outward square-root enclosure on request. General route: `boundary_distance` and `boundary_distance_with` combine `point_to_boundary` and `nearest_face` into a `BoundaryDistance` record with unsigned distance bounds, admitted witness, convergence and `BoundaryDirection` (`BoundaryNormal`, `UniqueNearest`, `WitnessOnly` or `Undefined`). | Fast box construction refuses non-box geometry. `BoundaryDirection` distinguishes an outward boundary normal at zero distance, a witness ray with a proved unique owner and separation margin, and a witness ray with tied contenders (such as outside a convex edge). It searches named faces only, so `distances()` refuses `solver.distance-coverage` when a live gyroid, sheet, lattice or mesh leaf has no named boundary face (regression: radius-2 sphere with period-4 level-1 sheet gyroid, boundary point (1,0,0)). |
| Q5 `distances().nearest_face`, `named_point`, `box_patch`, `face_patch` | Existing F8 ownership gives persistent origin, admitted trimmed witness, separation margin lower bound and contenders. `prepare_nearest` is called before exposing the distance service. Named point membership uses exact result trimming. `face_patch` and `face_patch_with` integrate `Fixed`, uniform `Traction`, area-normalized `TotalForce` and inward normal `Pressure` on any surviving named planar or axial-curved face (including holed plates, slotted walls, cylinders, cones, spheres, tori and split multi-piece faces), returning `FacePatch` with certified area, centroid, outward normal integral, resultant force and resultant moment enclosures. `box_patch` retains the exact rational shortcut on a single native axis box. | Ties and unsupported trims preserve F8 ambiguity or refusal. `face_patch` refuses absent names with `solver.face-name`, consumed or buried faces with `solver.patch-membership`, and non-axial curved trims outside the analytical integrator with `solver.patch-coverage`. |
| Q6 `material`, `Materials::Recipe` and Q3 pieces | Explicit homogeneous region, ordered axis slabs, or body-step material regions declared on `Recipe::materials` via `Materials::Recipe` (`Definition::from_recipe`). Regions intersect the part, never its nearest surface. Every piece carries region and material IDs; internal planar and curved interfaces retain adjacent identities in `MaterialAt::Interface`. | Gaps, overlaps, duplicate or empty names, and invalid axes or endpoints refuse `solver.material-partition` at compile or cell integration. Part boundary and unresolved point membership remain explicit. Constitutive parameter validation belongs to the solver. |

The fast box distance construction considers six faces for an interior point,
retaining every tied minimum, or clamps an exterior point to the closed box.
Thus squared length and displacement are exact rational constructions. Dividing
a nonzero displacement by the square root of its squared length defines its unit
direction; callers must preserve square-root enclosure error. This interface
returns that exact ray representation rather than an unproved floating normal.
A unique closest point at an exterior edge/corner still does not establish a
unique outward surface normal. `Query::direction` and F8's witness direction
retain their existing, distinct meanings.

The `_with` services require a preparation from the same `Definition` instance.
They reject another instance with `solver.definition-mismatch` before querying
geometry or stamping its revision. The definition's plan is private, so its
prepared service also retains the solver's named-boundary coverage check.

Cone and sphere integration checks whole-band trim coverage before omitting
an inactive band. Mixed full and partially surviving bands refuse
`solver.patch-coverage` when general quadrature is needed. Box material regions
must partition the complete certified bounds, even for a curved part. Gaps or
overlaps in that cover refuse `solver.material-partition`; the cell's certified
occupancy cannot be inferred from an incomplete material assignment.

## Integration proof and widths

For every cell, certified Inside pieces form a subset of its material, while
Inside plus Unknown pieces cover it. Outside pieces have zero material measure.
The pieces have disjoint interiors; shared facets and material interfaces have
zero volume. Exact box clipping intersects the solid and cell analytically.
Consequently summing piece volumes gives lower/upper material volume bounds,
and division by the positive exact cell volume gives an enclosure in `[0,1]`.
The grid partitions the complete certified compact bounds, so summing its
fractions times cell volume encloses the whole part volume. This is an enclosure
argument, not a point-sampling or corner-sign assumption.

For `xi_i=(x_i-cell_lo_i)/(cell_hi_i-cell_lo_i)`, return
`M_abc=integral xi_x^a xi_y^b xi_z^c dV`, with `0<=a,b,c<=2`, index `9a+3b+c`.
Each monomial is nonnegative on the cell. Its integral over a rectangular piece
is the product of three exact one-dimensional antiderivatives. Certain pieces
contribute that value to both ends; unknown pieces contribute `[0,value]`.
Summing by region encloses all 27 moments, including mixed degree-four and
-six terms required by Q1. `fraction.width()` times cell volume is the unresolved
volume; each moment has its own achieved `width()`. No small-cell clipping,
negative-value repair, discarded region, or `alpha * full_cell_stiffness` occurs.
An absent region query returns `solver.region-absent`, not an invented identity.

## Parallel work and limits

`grid` requires an active caller operation allowance. Default is a 64³ grid
(262,144 cells) over the certified compact bounds. It creates canonical jobs of
128 consecutive cells; x varies fastest, then y, then z. Each job has fresh
numeric caches and scratch. `OperationAllowance::reserve` charges the entire
allocation once; `split_off` transfers non-cloneable pieces to jobs. No worker
creates a replacement request allowance, refunds unused work, or borrows work
from another job. All jobs finish before the first canonical refusal is returned.
The default reservation is 256 units per cell, 67,108,864 total, plus 262,144
coordinator reduction units on success. This operation fits the 100M contract;
it refuses a default 64³ request under 1M or 20M. Those limits are unchanged.
A caller may explicitly request a smaller grid/reservation; insufficient job
allowance refuses rather than returning a partial enclosure.

Default workers come from the shared field pool, which respects
`available_parallelism()` and an enclosing request worker limit. Nonempty jobs
also cap the count. Explicit 1/2/4 caps the pool at that many workers when enough
jobs and processors exist. Work,
canonical cell order, all fractions, revisions, volume bounds and complete
Refusal equality (including attached undecided evidence) are compared directly,
not by Debug text. Worker count and wall time are separate execution metadata.
Finite tests are not universal scheduling evidence. Pool panics are
explicit refusals, outside the deterministic numeric contract. Caller deadlines
are checked before/after the bounded job phase; worker jobs are bounded by their
transferred allowances, not by a propagated wall-clock deadline.

## Cantilever demonstration

`solver_beam` uses the existing public 10×10×100 mm beam recipe and unchanged
example-private Q1 elasticity solver. F7 admits all 640 cells of the 4×4×40 grid
as exact full box intersections, verifies all 27 moments, identifies the explicit
material, integrates `beam/start` as Fixed and `beam/end` as a 100 mm² traction
patch. Traction `(0,-0.1,0)` gives `(0,-10,0)` N. The fixture supplies
`E=200000 N/mm², nu=0.3`; the adapter does not infer these from a material name.
All numerical displacements, residuals, reactions and stresses remain UNPROVED.
This is the previously frozen beam fixture, not the different beam proposed in
the broader design. It does not solve a sliver, plate, imported part, curved
boundary load, or heterogeneous material. Geometry-only sliver controls exercise
fractions 1/2, 1/16, 1/256 and 1/4096 and tangential shifts; they do not establish
conditioning/stabilization acceptance.

## Reproduction and evidence

Run from the repository root with ordinary Rust prerequisites:

```sh
cargo test -p catallax-field --test solver_queries
cargo run --release -p catallax-field --example solver_grid
cargo run --release -p catallax-field --example solver_beam
cargo clippy -p catallax-field --all-targets -- -D warnings
cargo fmt --all --check
git diff --check
```

`solver_grid` reads only public `parts/P1-bracket.json` and
`parts/P4-shelled-tray.json`, separately builds the B-rep and requests its volume
under the existing 100M limit, then checks the full independent B-rep enclosure
against the 64³ field volume interval at 1/2/4/default workers. Compile,
preparation, independent B-rep reference and each grid are explicitly separate
requests; their work is reported separately, not hidden inside a reset allowance.
The two-second target is not an acceptance gate; elapsed times are author-only
observations, not isolated performance benchmarks. Hosted integration and saved
answer preservation gates remain separate; no local protected replay was run.

## Executed first-cut results

Baseline: `39aa85a8513ca1b0283a5f28c188b74df6d08412`.
Raw public stdout is retained in [grid evidence](evidence/f7-grid.txt) and
[beam evidence](evidence/f7-beam.txt). The grid asserts enclosure of both ends of
the independent B-rep interval, and direct equality of complete ordered grid
records across actual 1/2/4/default worker runs. Default discovered four cores.

| Public part | Independent B-rep volume, approximately mm³ | Exact 64³ field volume interval, mm³ | Interval width / B-rep volume | Actual job work | Charged request work |
| --- | ---: | --- | ---: | ---: | ---: |
| P1 | 19943.937979737 | `[71425125/4096, 172245375/8192]` | 17.991793% | 1,872,688 | 67,371,011 |
| P4 | 4831.253636667 | `[263025/64, 42489225/8192]` | 22.290433% | 2,329,192 | 67,371,011 |

The charged work includes the fixed reservation, canonical reduction and
cancellation checkpoints. Unused reservation is not refunded. These intervals
are too wide to claim solver-quality integration. P1/P4 moment attribution was
also executed at every 4093rd cell: 12/59 nonempty cells respectively, checking
all 27 ordered intervals and equality of M000 with fraction times cell volume.
That is selected-cell execution evidence; the per-piece integration argument
above, not sampling, establishes the enclosure contract.

At the exterior probe `x=bounds.x_max+1`, y/z at the bounding-box centre, the P1
trimmed-distance and ownership probe completes in 9,808,982 work units under the
unchanged 100M (and 20M) allowance. Its boundary search converges in 32 cells to
`[427474432952006640711211374577123365817/42535295865117307932921825928971026432, 427474432952006640711211374577123365819/42535295865117307932921825928971026432] mm`
(around `sqrt(101) mm`, witness `(60, 20, 5)` on the shared edge between
`plate/side.right` and `plate/end`), and `nearest_face` reports the exact
two-contender tie (`field.ownership-ambiguous`, spread `1 / 2^127 mm`). P4's
probe uses 5,305,440 work units (down from 20,412,145); its boundary search
converges in 3 cells to `[1, 1] mm`, and its owning-face query proves
`blank/x_max` with exact owner distance `[1, 1] mm`, a witness ray along
negative x and separation margin lower bound `1/2 mm`.

The F7 beam adapter matches the published frozen 4×4×40 tip
`-0.01929718439662 mm`, 166 iterations and residual `7.160821669988e-10`.
Its Euler–Bernoulli difference is 3.514078%; this grid does not meet M20's 2%
beam mark. Geometry/integration preparation uses 22,400 work units. Stress and
reaction rows remain UNPROVED; no new reaction determinism claim is made.

Thirteen focused Rust tests pass in `solver_queries`, covering the analytic moment
tensor, geometry-only slivers, unknown sphere moment enclosure, named box and
general planar/curved face patches (`plate/start` with bore and bolt holes,
`flange/side.top` with rounded slot, `boss/side.wall`, `slot-tool-0/side.wall`,
and P4's four-piece `rounding/side.wall` and trimmed `blank/side.x_max`),
general `boundary_distance` and P1/P4 distance probes across 1/2/4 workers,
the unnamed sheet-boundary distance refusal, slab and curved recipe-defined
material regions (`Materials::Recipe` with a cylindrical copper core inside a
steel shell), mixed partial-cap refusal, narrow material gaps and overlaps in a
sphere, cross-definition preparation rejection, and direct semantic comparisons
at 1/2/3/4/default workers. The curved material grid has 512 cells in four jobs
and asserts the actual worker count for its 1/2/4-worker comparisons. Both
caller-reservation exhaustion (unchanged 1M and 20M limits) and child-job
exhaustion retain equal refusals and work. A 300-region control proves
material partition scans cannot escape job accounting. All-target Clippy with
warnings denied, whole-workspace formatting check and diff check pass. Cargo
emits its existing duplicate-bin manifest warning. No hosted gate, protected
saved-answer replay, C++ oracle, general sliver elasticity solve or M20 acceptance
is claimed by this evidence.
