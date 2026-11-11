# M10 cylinder and sphere sections

The second slice cuts one `RationalPatch` by a cylinder or a sphere with the
same complete search, exact ends, typed refusals and candidate tube as
[the plane slice](M10-PLANE-SECTIONS.md). Only the cut polynomial changes.
`SectionSurface` names the admitted cutting surface: plane, cylinder or sphere.
`RationalPatch::sections` takes one and returns `PatchSectionBranch` values,
retained as `Curve::PatchSection`. Cones and tori are later slices.

## The cut polynomial

The surface's implicit equation is completed to a homogeneous form in the
patch's `X, Y, Z, W` and substituted exactly. A plane stays linear. A cylinder
or sphere is a quadratic form, so a patch of bidegree `(m, n)` gives a cut
polynomial of bidegree `(2m, 2n)`. Its Bernstein coefficients over a span, or
over a tube box, come from the exact Bernstein product of the homogeneous
control nets: the coefficient at `(i, j)` sums
`C(m,i1) C(m,i2) C(n,j1) C(n,j2) B(h1, h2) / (C(2m,i) C(2n,j))`, where `B` is
the symmetric bilinear form and `i1 + i2 = i`, `j1 + j2 = j`. Every rational
operation is reserved against the caller's operation allowance.

On the surface the gradient of the cut polynomial is `W^d` times the surface
gradient composed with the patch derivatives. Weights are positive, so its
signs are the surface's. The cell tests, exact edge roots, cell joins,
traversal, membership, order, tangents and pair orientation therefore apply
unchanged. Orientation follows the cutting surface's outward normal cross the
patch normal, and reverses when the input surfaces swap. Reducible spans keep
the exact rational spline reductions of the first slice.

## Candidate tubes

The tube checker works on the cut polynomial, so cylinder and sphere tubes use
the same walls, transverse slope bound, end caps and speed bound. One rule
changed. The first checker required every arc cell to have been classified
along one parameter. A quadric cut can be proved monotone in `v` in one cell
and in `u` in the next. The checker now accepts an open single-span branch
whose ends lie on opposite borders of either parameter, and whose cell ends
run monotonically in it. The tube still proves its own transverse uniqueness
over every candidate box, so the cell classification never enters the
theorem. The plane slice's decisions run in the same order, and its recorded
work units are unchanged.

The checker refuses a candidate whose end control is not exactly on the
border. The tight SMLib sphere trace ends 4.3e-19 off the `u = 0` border, and
its refusal is a regression case. The test caller moves end controls within
1e-12 of a border onto it. The checker treats the result as one more
untrusted candidate.

## Verification

```sh
cargo test -p catallax-intersect --test patch_quadric -- --nocapture
```

The loft-like patch of the first slice is cut by the cylinder
`x² + z² = 25/4` and by the sphere `x² + y² + z² = 25/4`. Both cuts end
exactly at `(0, 0, 5/2)` and `(2, 0, 3/2)`. The hill is cut by a radius-1
cylinder along `x` and a radius-1 sphere above the hilltop. Each touches the
hilltop exactly, sits one millionth below it, one millionth above it, or has
its lowest point at 0.3.

The pinned SMLib returns no curve for either surface one millionth below the
hilltop. The complete search finds a closed loop in both cases. At the exact
hilltop it refuses with `intersect.patch-cylinder-tangent` or
`intersect.patch-sphere-tangent` and names `(5/2, 5/2, 225/256)` at parameters
`(3/2, 3/2)`. SMLib returns nothing there either. Fixture provenance is in the
[oracle fixture notes](../crates/catallax-intersect/tests/fixtures/patch-quadric/README.md).
Raw evidence retained privately.

The measured search operations, including patch admission, are below. The
allowance remains 20,000,000 units for every row.

| Case | Result | Operation units |
| --- | --- | ---: |
| Cylinder, loft | 1 open branch, 2 arc and 1 empty cells, depth 2 | 106,943 |
| Cylinder, hill at 0.3 | 1 closed loop, 10 arc and 5 empty cells, depth 2 | 1,102,170 |
| Cylinder, 1e-6 below the hilltop | 1 closed loop, 12 arc and 80 empty cells, depth 22 | 928,053 |
| Cylinder, exact hilltop | Typed tangent refusal naming the touch point | 209,226 |
| Cylinder, 1e-6 above the hilltop | Proved empty, 12 empty cells | 334,913 |
| Sphere, loft | 1 open branch, 2 arc and 1 empty cells, depth 2 | 148,707 |
| Sphere, hill at 0.3 | 1 closed loop, 12 arc and 10 empty cells, depth 3 | 963,754 |
| Sphere, 1e-6 below the hilltop | 1 closed loop, 8 arc and 84 empty cells, depth 22 | 889,658 |
| Sphere, exact hilltop | Typed tangent refusal naming the touch point | 266,895 |
| Sphere, 1e-6 above the hilltop | Proved empty, 12 empty cells | 438,706 |

Candidate checks on the loft cuts, including admission and search:

| Candidate | Result | Operation units |
| --- | --- | ---: |
| Cylinder, tight (16 pieces) | U < 91.871 nm, rho = 1/65536 | 1,147,992 |
| Cylinder, default (7 pieces) | Accepted under the 1 µm cap only, U < 735.0 nm, rho = 1/8192 | 661,183 |
| Cylinder, tight shifted by 1e-3 | Typed cap refusal | 107,981 |
| Cylinder, tight stopping short | Typed endpoint refusal | 106,955 |
| Sphere, tight as traced (17 pieces) | Typed endpoint refusal, end 4.3e-19 off the border | 148,719 |
| Sphere, tight with ends on the border | U < 91.871 nm, rho = 1/65536 | 3,759,489 |
| Sphere, default (8 pieces) | Typed cap refusal | 3,497,941 |
| Sphere, tight shifted by 1e-3 | Typed cap refusal | 149,745 |
| Sphere, tight stopping short | Typed endpoint refusal | 148,719 |

The 1-, 2- and 4-worker checks agree on exact results, endpoint enclosures,
refusal details, candidate bounds and operation units. With their endpoint
enclosures and tight candidates the loft rows cost 1,148,088 units for the
cylinder and 3,759,585 for the sphere. The hill rows cost 1,103,386 and
964,937, and the tiny loops 929,364 and 890,969.

Pair queries also run for every case and both input orders. The largest,
the sphere tiny loop with the patch first, uses 1,553,784 units. Its exact
check that each returned tangent is perpendicular to the surface gradient
runs outside the query's allowance. At algebraic loop witnesses the default
decider bounds that product below 1e-40 instead of proving it zero. A Boolean
that consumes these branches will meet the same cost.

Atlas joins, third-surface meetings and the first Boolean consumption are in
[M10 atlas sections and meetings](M10-ATLAS-SECTIONS.md). There are no
spline-against-spline intersections or STEP writer for these branches yet.
Closed and multi-span branches keep their exact definitions; the candidate
checker refuses them.
