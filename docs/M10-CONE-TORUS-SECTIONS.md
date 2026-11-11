# M10 cone and torus sections

The third slice extends [cylinder and sphere sections](M10-QUADRIC-SECTIONS.md)
to a cone's retained nappe and a ring torus. `SectionSurface::Cone` and
`SectionSurface::Torus` use the same complete cell search, exact endpoints,
component traversal, parameter queries and candidate tube certificates.
The search runs without a candidate. Every accepted branch retains its patch,
cutting surface and complete cell cover as `Curve::PatchSection`.

## Cut equations and completeness

A cone uses the existing homogeneous quadratic substitution. Its polynomial
contains both nappes. Before either a rational reduction or a cell search,
the homogeneous height controls must prove the whole span strictly inside
the retained nappe or strictly outside it. An outside span is empty. An
unresolved height net, including an apex or a nappe boundary, refuses with
`intersect.patch-cone-nappe`. This conservative rule does not clip a span
across the apex plane.

For a torus, let `D = (X,Y,Z) - center W`, `Q = D dot D`,
`H = D dot axis`, and `N = axis dot axis`. Its homogeneous cut is

```text
A = Q + (major² - minor²) W²
B = Q - H²/N
F = A² - 4 major² B W²
```

The existing quadratic Bernstein product builds `A`, `B` and `W²`.
Two further tensor Bernstein products build the quartic. A patch of bidegree
`(m,n)` gives a torus net of bidegree `(4m,4n)`. A polynomial identity test
checks both this net and the span equation against the homogeneous global
implicit equation, with nonconstant weights, a translated centre and an
oblique, nonunit axis. Only ring tori with `0 < minor < major` are admitted.

Positive patch weights preserve the sign of each homogeneous cut. The search
therefore uses the earlier slices' empty-cell and unique-arc proofs without
changing them. Exact boundary roots join every accepted cell into complete
open components or closed loops. A candidate tube can represent an already
proved component, but cannot establish that no component was missed.

The tangent-touch prototype remains parked. This slice accepts regular cuts
only. An exact stationary contact or a repeated root in a rational reduction
refuses with `intersect.patch-{cone,torus}-tangent`. An unresolved cell search
refuses with `intersect.patch-{cone,torus}-nonregular`. The boundary-contact
case below reaches that named refusal at the unchanged depth allowance.
No tolerance, depth, cell or operation allowance changes.

## Verification

```sh
cargo test -p catallax-intersect --test patch_cone_torus -- --nocapture
cargo test -p catallax-shapes quartic_net_and_span_equation
```

The open patch is `z = 2 + (x² - 9/25)/8`, with `x` between `-3/5` and
`3/5`, and `y` between `7/10` and `3/2`. A cone of slope `1/2` about `z`
and a torus centred at `(0,0,2)` with radii `3/4` and `1/4` each cut one
open component. Both have exact endpoints `(-3/5,4/5,2)` and `(3/5,4/5,2)`.

The dome is `z = 1 - x² - y²`, with both coordinates between `-1/4` and
`1/4`. For an axis through the dome, the cone has slope `1/10` and apex
at the origin. The torus has centre `(0,0,1)` and radii `1/8` and `1/32`.
The cone cuts one loop; the torus cuts two nested loops.

The near-tangent cone has apex `(-1,0,low)`, axis `(1,0,1)` and slope `1`.
The near-tangent torus has centre `(1,0,low + 1/4)`, axis `z` and radii
`1` and `1/4`. Each cuts a loop at `low = 99/100`, a tiny loop at
`1 - 1e-6`, touches exactly at `1`, and misses at `1 + 1e-6`.
The exact contact refusal names `(0,0,1)` at patch parameters `(1/2,1/2)`.

Search measurements include patch admission and endpoint equality checks for
the open cuts. The allowance is 20,000,000 operation units for every row.

| Case | Result | Operation units |
| --- | --- | ---: |
| Cone, open cut | 1 open branch, 1 arc cell, depth 0 | 5,458 |
| Cone, axis through dome | 1 loop, 8 arc and 12 empty cells, depth 6 | 296,916 |
| Cone, low = 99/100 | 1 loop, 8 arc and 12 empty cells, depth 6 | 73,982 |
| Cone, 1e-6 below contact | 1 loop, 10 arc and 62 empty cells, depth 20 | 315,159 |
| Cone, exact contact | Typed tangent refusal | 15,782 |
| Cone, 1e-6 above contact | Proved empty, 4 empty cells | 18,237 |
| Torus, open cut | 1 open branch, 8 arc and 6 empty cells, depth 6 | 1,396,177 |
| Torus, axis through dome | 2 loops, 48 arc and 40 empty cells, depth 10 | 1,359,955 |
| Torus, low = 99/100 | 1 loop, 8 arc and 16 empty cells, depth 7 | 1,813,688 |
| Torus, 1e-6 below contact | 1 loop, 12 arc and 64 empty cells, depth 20 | 1,720,555 |
| Torus, exact contact | Typed tangent refusal | 50,595 |
| Torus, 1e-6 above contact | Proved empty, 4 empty cells | 78,090 |

Sixteen cubic Hermite pieces, computed in the test and retained as untrusted
binary64 candidates, represent each open cut. Each accepted tube has
`rho = 1/8192` and `U < 97.657 nm`, below the unchanged 100 nm target.
Admission, complete search and candidate checking cost 91,557 units for the
cone and 2,174,019 for the torus. Moving the candidate by `1/100` in its
transverse parameter gives a typed cap refusal. Stopping short gives a typed
endpoint refusal. A one-unit operation budget yields a typed budget refusal
with no partial cut.

Additional cases exclude the opposite cone nappe, refuse a span that reaches
the apex plane, and refuse flush cone-generator and torus-circle contacts.
The flush and boundary cone cases cost 1,917 and 1,914 units. The flush torus
circle costs 9,555 units. The isolated torus boundary contact costs 28,693
units and returns `intersect.patch-torus-nonregular`.

Every regular-search case and every flush or boundary-contact case also runs
after translation, an oblique rational rotation and translation, and scales
of `1/10` and `10` followed by that pose. The regular cuts keep the same cell
counts, closure and parameter-endpoint enclosures. Open endpoints agree
exactly after undoing each copy's transform. Contact copies keep the same
named refusal, and isolated tangent refusals name the transformed touch point.
The whole matrix agrees at 1, 2 and 4 workers, including each case's work and
refusal detail. Copies can cost different work; no wall-clock limit gates
this capability. The largest search-and-endpoint row is 1,825,259 units.

Pair queries run each regular case in both input orders. They check branch
membership and exclusion, open order, closed cyclic order, nonzero tangents,
and orientation under swapping the surfaces. Closed-branch candidate checks
refuse with the family-specific domain code. The largest pair query is the
torus tiny loop with the patch first, at 2,410,572 units. Independent tangent
checks enclose the gradient product below `1e-40`; those checks run outside
the query allowance.

The earlier slices' limits remain. Atlas joins, third-surface queries,
spline-against-spline cuts, Boolean consumption and STEP writing for retained
branches remain open. Closed and multi-span branches retain exact definitions;
the candidate tube checker still requires one open single-span graph branch.
Raw evidence retained privately.
