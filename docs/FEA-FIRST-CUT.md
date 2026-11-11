# Unproved rectangular-beam experiment

This example solves three-dimensional, small-strain, homogeneous isotropic
linear elasticity on a plain Cartesian grid. **Every displacement, stress,
reaction and residual is an unproved numerical estimate.** It is a first-cut
experiment, not production FEA, a stress certificate, an export qualification,
or completion of M16/M20. It has no route into proven answers or export.

The draft [FEA design in PR385](https://github.com/BohnBawerick/catallax/pull/385)
at `0f6ae66a2a1d5481995f390fca1f25383bc54690` is context only. This example does
not adopt its unresolved owner decisions, cut-cell method, material-region
interface or proposed benchmark geometry. In particular, this experiment uses
a 100 by 10 by 10 mm beam, rather than the draft's 100 by 5 by 5 mm beam.

## Reproduce

From the repository root:

```sh
cargo test --release -p catallax-field --example beam_first_cut
cargo clippy -p catallax-field --example beam_first_cut -- -D warnings
cargo run --release -p catallax-field --example beam_first_cut -- 2
cargo run --release -p catallax-field --example beam_first_cut -- 4
cargo run --release -p catallax-field --example beam_first_cut -- 8
```

The optional second argument is a recipe filename. The public input is
[`beam/recipe.json`](../crates/catallax-field/examples/beam/recipe.json).
Only this exact parsed recipe is admitted: unsupported geometry, changed
placements/dimensions or naming refuse as `fea.unsupported-geometry`. The only
supported boundary conditions are the named clamp and load below; there is no
general constraint input. Extra command arguments fail. Grid choices other
than 2, 4 or 8 cross-section cells refuse as `fea.grid`.

## Frozen problem and method

* Geometry: x and y in [0,10] mm, z in [0,100] mm, a single recipe block.
* Material: E = 200000 N/mm², Poisson ratio = 0.3; no density/body force.
* Support: all three displacement components zero on `beam/start`, z = 0.
* Load: uniform traction (0,-0.1,0) N/mm² on `beam/end`, z = 100,
  resultant (0,-10,0) N, with no imposed end moment. Other faces are free.
* Grids: [2,2,20], [4,4,40], [8,8,80], with cubic cells of 5, 2.5, 1.25 mm.
* Elements: continuous trilinear Q1 displacement, three components per node;
  full 2 by 2 by 2 Gauss integration of `Bᵀ C B`, no reduced integration,
  fictitious stiffness, beam kinematics or empirical bending correction.
* Elasticity: engineering strains ordered xx, yy, zz, xy, yz, xz. Normal
  constitutive block has lambda on every entry and 2 mu on its diagonal;
  shear entries are mu. Lambda and mu are the ordinary isotropic Lamé values.
* Fixed degrees of freedom are eliminated symmetrically. Each loaded surface
  cell contributes one quarter of its traction resultant to each face node.
* Fixed-order sparse assembly and Jacobi-preconditioned CG; relative residual
  target 1e-9, maximum 10000 iterations, maximum 5120 cells. The original
  residual `f-Ku` is recomputed before accepting convergence. Exhaustion,
  nonpositive curvature or nonfinite quantities return `fea.solver`.

The immutable recipe whitelist establishes the aligned full-cell domain; finite
point samples are **not** a proof of arbitrary cell occupancy. The example
compiles the real recipe to a field, checks every cell centre using
`proven_checked`, and checks outside halo probes. It selects both origins with
`Plan::select` and checks every end patch centre with `face_query_checked`.
Unresolved or unexpected classifications refuse. These queries share one
unchanged default Decider policy; no work allowance is reset per cell. No fast
preview signs feed assembly. No mesh or C++ oracle is involved.

## Observables and mathematical references

The tip observable is the area-mean y displacement on the loaded face, using
the same consistent weights as the traction. Euler–Bernoulli gives
`I = b h³/12 = 833.333333… mm⁴` and
`delta_y = F L³/(3 E I) = -0.02 mm`.
This formula is a comparison only; it is not used in assembly or solution.

Stress is recovered from `C B u` at cell centres, without smoothing or
extrapolation. The report samples one line along z at x = 5 + dx/2 and
y = 10 - dy/2. It compares von Mises stress with the absolute bending stress
`|F (L-z) (y-5)/I|` at those same coordinates. The von Mises calculation uses all
six stress components. Beam bending alone omits transverse shear and 3D end
constraint effects; agreement is not expected at either end. Sampling points
approach the surface as the grid refines; they are not boundary peak stresses.
The clamp can have stress singularities; no maximum-stress claim is made.

Support reactions use the original element stiffness before elimination. The
expected net reaction is (0,10,0) N and expected moment about the root centre is
(-1000,0,0) N mm. Equilibrium and a small algebraic residual do not certify the
continuum discretization error.

## Validation and limits

The example tests affine uniaxial energy against an independently calculable
value, element symmetry and rigid translations, load reversal and scaling,
grid refinement, force/moment balance, exhaustion and invalid inputs. Independent
1/2/4-worker solves compare displacement, residual, iterations and stress bits;
the solver itself remains serial with fixed reduction order. This does not
claim portability of binary64 estimates across different processor platforms.

No new shared API, dependency, production module, saved answer, tolerance or
existing allowance is changed. Other cloud lanes have no API dependency on this
example. Private inputs and saved-answer inventories are neither read nor run.
Existing hosted CI remains responsible for integration and preservation gates;
focused Rust validation is not an overall PASS claim.

Unsupported: arbitrary recipes, cut cells, rotated grids, imported geometry,
other supports, material regions, contact, large deformation, plasticity,
buckling, dynamics, nearly incompressible materials, certified stress bounds,
and a UI stress picture. Grid refinement does not establish any of these.

## Recorded run

Baseline: `371774eaf3ffb8321ad54311c7596a25c8a0cce3`. The following are
**unproved** release-build estimates on the declared public recipe. The
Euler–Bernoulli reference is -0.020000000000 mm for every row.

| Grid | Tip mean y (mm) | Absolute relative difference | CG iterations | Actual relative residual |
| --- | ---: | ---: | ---: | ---: |
| 2 × 2 × 20 | -0.01751578866283 | 12.421057% | 79 | 7.32343e-10 |
| 4 × 4 × 40 | -0.01929718439662 | 3.514078% | 166 | 7.16082e-10 |
| 8 × 8 × 80 | -0.01982293562727 | 0.885322% | 342 | 9.02912e-10 |

All three underestimate the displacement magnitude. The finest experimental
configuration is within 2% of its beam formula; this is **not** a pass of the
proposed PR385 geometry or overall M20. The refinement test requires decreasing
deflection error and the finest error below 2%; the roadmap mark is not relaxed.

Representative unproved stress comparisons (N/mm²) near midspan:

| Grid | z (mm) | y (mm) | von Mises | Beam bending magnitude | Relative difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 × 2 × 20 | 47.5 | 7.5 | 1.402082896623 | 1.575000000000 | 10.978864% |
| 4 × 4 × 40 | 48.75 | 8.75 | 2.235022194264 | 2.306250000000 | 3.088469% |
| 8 × 8 × 80 | 49.375 | 9.375 | 2.636668825480 | 2.657812500000 | 0.795529% |

The maximum sampled relative stress difference over 20 ≤ z ≤ 80 mm is
11.358413%, 3.148221%, and 0.808364%, respectively. This interior interval is
only a reporting summary, not a stress acceptance threshold. Every cell-centre
row, including both end regions, is emitted by the example. At the finest first
cell, z = 0.625 mm, von Mises is 4.055653705783 versus beam bending
5.217187500000 N/mm²: the end discrepancy remains visible.

The finest net y reaction is 10.000000000907 N and root x moment is
-1000.000000008875 N mm. Other reaction force components have magnitude below
3.1e-9 N and other moment components below 2.6e-8 N mm.

Focused release tests (four), all-target field Clippy with warnings denied,
workspace formatting and whitespace checks passed. No full workspace test,
private replay, C++ check, imported benchmark or hosted integration result is
claimed here. Hosted CI must validate the final PR head.
