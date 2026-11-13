# Agreed M20 cantilever

`m20_beam` runs the first rectangular cantilever in the agreed
[FEA design](FEA-DESIGN.md). Displacement, stress, reactions and residuals are
unproved numerical estimates. Geometry certificates do not certify elasticity.
The beam, plate, sliver and imported-part acceptance lines remain separate.

## Reproduce

```sh
cargo test --release -p catallax-field --example m20_beam --example solver_beam --example beam_first_cut
cargo run --release -p catallax-field --example m20_beam -- m20-output
cargo clippy -p catallax-field --example m20_beam --example solver_beam --example beam_first_cut -- -D warnings
cargo fmt --all --check
```

The output directory contains one directory per grid. Each has `beam.svg`,
`profiles.csv` and `cell-stresses.csv`. The SVG shows cell-centre von Mises
stress on a declared longitudinal slice and the centreline displacement curve
beside Euler-Bernoulli. Its thickness is magnified and labelled. It does not
smooth stresses or substitute beam theory for the field solution.

## Geometry and solver policy

The public input is [beam/agreed.json](../crates/catallax-field/examples/beam/agreed.json).
It is a single box with world bounds `[0,100] x [-2.5,2.5] x [-2.5,2.5] mm`.
The recipe independently builds a B-rep with exact volume `2500 mm³`.
M16 integration admits each assembled cell only with exact full material
fraction, a connected box and 27 exact material moments. The accumulated volume
must equal `2500 mm³`. There are no preview signs or sampled occupancy decisions.
Both representations and the picture describe this same rectangular recipe.
This run does not establish export agreement for curved or imported parts.

The fixed face is `beam/side.x_min`. The loaded face is `beam/side.x_max`.
M16 resolves both names at their world positions and checks each area is
`25 mm²`. The explicit traction is `(0,0,-0.04) N/mm²`, with resultant
`(0,0,-1) N`. Material assignment is homogeneous, with explicit
`E=200000 N/mm²` and `nu=0.3`. Other faces are traction-free.

The existing example-private solver uses local axes `(world y, world z, world x)`.
This cyclic permutation is a proper rotation. Output maps vectors, moments and
all six stress tensor components back to world axes. The fixed world x face
eliminates all displacement components. The loaded face uses consistent Q1
surface weights. Integration is full `2 x 2 x 2` Gauss integration of `Bᵀ C B`.
There is no beam element, reduced integration or fitted stiffness correction.

The frozen world grids are `[20,2,2]`, `[40,4,4]` and `[80,8,8]`.
The solver keeps its existing limits of 5120 cells, 80 cells per axis,
10000 iterations and verified relative residual at most `1e-9`.
Jacobi-preconditioned CG and its fixed assembly/reduction order are unchanged.
Each complete geometry request retains the existing 100M operation allowance.
The report separates charged geometry work, iterations, matrix entries,
matrix-vector products and their scalar products. These solver counts do not
claim to be a certified numeric operation allowance or a timing measurement.

Each grid runs complete independent requests concurrently with 1, 2 and 4
workers. Direct equality checks include every displacement unknown, recovered
stress row, reaction, residual, iteration count, solver count and charged
geometry work. The elasticity solver itself is serial. This tests concurrent
request determinism, not parallel assembly or cross-processor bit equality.

## Stress recovery and comparison

`profiles.csv` includes world x stations every `1.25 mm`, on both outer-fibre
paths `(x,0,+/-2.5)` and the centreline `(x,0,0)`. Each incident material cell
provides its own `C sym(grad(u))` tensor and full-tensor von Mises value.
At cell interfaces the CSV retains every trace; stdout gives their range at
`x=25,50,75 mm`. The corresponding outer-fibre beam references are
`3.6,2.4,1.2 N/mm²`. Centreline bending reference is zero and appears separately.
End-region traces remain in the CSV and are outside these span comparisons.
There is no separate stress-percentage acceptance claim.

The profile includes signed and absolute displacement and stress differences.
`cell-stresses.csv` contains every cell-centre tensor, without smoothing.
Tensor order is `xx,yy,zz,xy,yz,xz`, in `N/mm²`. Displacements use mm.
An affine patch test checks stress recovery and displacement interpolation at
interior points and all incident-cell traces on internal and exterior facets.

The deflection observable is the consistently integrated area-mean tip
displacement. Its independent Euler-Bernoulli reference is `-0.032 mm`.
The comparison returns the example-private typed refusal
`fea.deflection-accuracy` whenever `abs(tip / -0.032 - 1) > 0.02`, including the
attained displacement, reference and error. A refusal preserves the numerical
field and picture, labelled unproved. It does not qualify a missed mark as a pass.

## Remaining M20 work

The plate with a hole, aggregate-constrained sliver solve and complete imported
solve remain open. M16's exact sliver measures do not demonstrate sliver
elasticity. The plate stress rule and imported configuration follow the
sequencing in D5 and D6. No tolerance, work limit or saved replay answer changes.
Raw evidence retained privately.
