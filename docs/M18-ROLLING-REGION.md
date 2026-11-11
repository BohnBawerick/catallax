# M18 finite lower rolling-ball region

**Partial implementation; authored P5 and section 4.17 smooth blend remain open.**
This continuation implements an exact, independently queryable lower-envelope
region. It does not change recipe fillet admission. The new region's set is
explicit below: equivalence to the complete trimmed recipe fillet is still an
obligation, not an assumption or a sampled result.

The dependency is the finite-clearance groundwork at
`cdd11fddbfd99855224f7b0b655fd08961c7a1e0` (PR #409). Fresh main was
`1a334a42e6f70b04cdd62bdfd920d6f8441cd9b0`. The continuation merges these at
`6845e7827b850e2d4ad9614554c73f180b936572`; the groundwork branch is unchanged.

## Exact set and coverage

`catallax_field::rolling_region::Region` consumes a validated
`RollingEnvelopeDefinition`, not a finished B-rep or a fitted carrier.
In its orthonormal frame, let `R` be body radius, `r` ball radius, `m` mouth
radius and `o` mouth-axis offset; `C=R+r`. The finite centre curve is

```
cos(theta0) = (C*C + o*o - m*m)/(2*C*o)
t0 = sin(theta0)/(1+cos(theta0))
t = u*t0,                        -1 <= u <= 1
cx = C*(1-t*t)/(1+t*t)
cz = 2*C*t/(1+t*t)
D = sqrt((cx-o)^2 + cz^2)
H = sqrt(r*r - (m-D)^2)
```

All quantities are exact expressions; all query arithmetic is outward rounded.
The definition constructor checks the rim normal-cone regime. The region also
requires an acute, nonempty mouth sector. It is clipped to
`R <= sqrt(x*x+z*z) <= C`, `0 <= y <= r`, `x >= 0`, and
`|z|*cos(theta0) <= x*sin(theta0)`.
For **every** centre parameter `u` whose planar distance `d` from `(x,z)` is at
most `r`, it additionally requires

```
y <= H(u) - sqrt(r*r-d*d).
```

This defines the material below the lower envelope, rather than the complement
of the balls on both sides. In particular, a point above a ball centre is not
admitted merely because it lies outside that ball. All lengths are millimetres.

`classify` proves strict membership in this explicit set. A midpoint with
`d^2<r^2` and either `y>H` or three-dimensional distance squared below `r^2`
is an exact exclusion witness. Inclusion requires a **whole interval** cover of
`[-1,1]`: each panel proves either `d^2>r^2`, or `y<H` and three-dimensional
distance squared above `r^2`, throughout that panel. Subdivision has no gaps;
shared endpoints are included. Sampling cannot establish inclusion.

The invariant `0 <= m-D <= m-(o-C)` tightens the interval evaluation, using only
bounds justified by the validated finite parameter domain. Arithmetic uses
96 bits, at most 4096 evaluated panels and depth 64. These are new local refusal
limits, not increases to any existing allowance or changes to geometric
accuracy. Unseparated boundaries, clipping-wall incidence and precision/resource
limits return typed refusals. There is no approximate success fallback.

`station` returns the exact centre, cylinder contact and mouth contact for any
rational normalized half-angle in `[-1,1]`, including both plane transitions.
No exact distance, unique direction, named recipe face, mass, print or export
certificate is implied by `Region::classify` or `station`.

## Executed controls

The controls use the public P5 dimensions (`R=12`, `m=4`, `o=17`, `r=2`, plane
height 5), with radius `5/2` contact controls, translation `(7,11,-13)` and a
90-degree rotation about z. Exact algebraic checks establish contact incidence
with both supports, radius equality, and both transitions' agreement with
`RollingEnvelopeDefinition::corners`.

A full-torus shortcut gives the wrong material at `(12.5,5.5,0)`:
its squared distance to the ordinary torus centre `(14,7,0)` is `4.5>4`, while
its squared distance to the native middle ball centre `(14,5+sqrt(3),0)` is
`5.5-sqrt(3)<4`. The new region correctly excludes it. `(12.2,5.2,0)` is proved
inside by complete interval coverage. The middle meridian has lower envelope
height `5+sqrt(3)-sqrt(7)/2` at `x=12.5`; controls at heights `5.4` and `5.42`
prove opposite sides. The exact boundary remains an explicit refusal.

These are analytic controls of the explicit native set, not a global
comparison against the B-rep. The radius `7/2` normal-cone exit still refuses.

## Production parallel service

`classify_batch` is an actual production method on `Region`, not a parallel test
wrapper around a serial API. Its default uses all available cores, capped by
input count. Fixed per-input allowances are transferred from a single
caller-owned `OperationAllowance`; unused units are not refunded. It inherits
the caller's numeric policy and absolute deadline (the only numbers API change
is a read-only `Decider::deadline` getter). Results and errors have input order.
Spawn failure/panic or caller cancellation returns a refusal without partial
success. Normal deterministic comparisons exclude wall-clock cancellation.

The seven-point matrix used actual 1, 2, 4 and default workers (default 4).
Four-worker partitions held 2/2/2/1 records. Full equality covers membership,
complete `Refusal` values including attached evidence, actual work and exhaustion.
At a 700000-unit reservation the point work was
`[184,156,156,6,156,2,2]` (662 total); caller accounting was 700002 including
its two cancellation checkpoints. Reservations of 7 and 0 produced identical
per-point exhausted refusals at all tested counts, with respective actual work
`[1,1,1,1,1,1,1]` and `[0,0,0,0,0,0,0]`. There is no worker-dependent budget
multiplication. Finite tests do not establish every possible thread count.

This does **not** parallelize recipe compilation. That remains a separate limit.

## Remaining blocker to the requested P5 exit

The explicit lower-envelope set has not been proved equivalent to the complete
recipe's finite contact region. Before recipe admission, the compiler still
needs all of these together:

- discovery of each finite mouth from the recipe's Boolean material, with
  exhaustive incidence of every selected support and every intervening feature;
- proof that the four regions, intervening torus sectors, mouth arcs and
  cylinder/plane rails cover exactly the required material, with no missing or
  overlapping transition and no interaction with another trimmed feature;
- exact boundary/name/ray integration for these regions. The field's current
  `Face` registry uses analytic support polynomials; substituting the fitted
  carrier, an arbitrary polynomial, or the complete torus is not valid.

The existing `fillet.unsupported-overflow` therefore remains for authored P5.
This is substantive native region work, **not** completion of the assigned P5
exit. Smooth blend is not implemented after that unresolved first exit. It
still requires true operand distances, not existing Boolean residuals. Wider
sweeps/lofts and complete M18 acceptance are also open.

## Reproduction and evidence limits

```sh
cargo test -p catallax-field --test rolling_region -- --nocapture
cargo test -p catallax-field --test m18_clearance --test fillets
cargo clippy -p catallax-field --all-targets -- -D warnings
cargo fmt --all --check
git diff --check
```

Six new tests and seven existing focused fillet/clearance tests passed. Strict
all-target field Clippy, formatting and diff checks passed. The
retained [worker records](evidence/m18-rolling-matrix.txt) contain the final
production batch outcomes. Raw local test/build logs are retained separately.
No C++ oracle, protected corpus, tolerance/allowance change, new recipe admission,
CI launch or hosted-gate inspection is included. Hosted integration is unverified.
