# M18 finite round-body clearance groundwork

This increment does **not** close the authored P5 field refusal and does not
implement the section 4.17 smooth blend. It removes a smaller, independently
provable obstruction found while investigating P5. The requested order remains:
first native P5 rim-contact continuation, then the separate smooth-blend feature.
No M18 completion or two-kernel parity is claimed.

Baseline: `a35b3f34d5b4031696d07dadfce32c123c41a165`.

## Change and proof

The old field fillet check rejects a candidate whenever its complete supporting
torus intersects another support and existing carrier/hull separation cannot
exclude it. For a finite countersink cone or cap, that can reject a valid fillet:
its infinite support crosses the torus but its finite source body misses the
whole material-change region.

The new fallback applies only after that existing intersection check would
refuse. It requires a concave circular plane/cylinder boss-foot fillet whose
torus major radius is exactly `R+r`, where `R` is the selected ring radius and
`r` is the rolling radius. Its short meridional arc and both contact rails bound
a material-addition region contained in the closed radial band `[R,R+r]`.
The other face must belong to a native round primitive (cylinder or capped
frustum) with a parallel axis. Every point of that primitive is within distance
`m=max(base_radius,top_radius)` of its axis. Its axis is at radial offset `d`
from the fillet axis. The complete finite primitive, hence every trimmed face it
owns, is strictly separated from the band if either:

- `d > R+r+m`; or
- `R > m` and `d < R-m`.

The implementation compares exact squared expressions, without rounded lengths,
using the existing rigid placements and exact sign decider. Equality remains
unresolved by this certificate. Other primitive kinds, nonparallel axes,
convex rings, hole-floor rings and actual intrusions retain the old path/refusal.
The proof is independent of feature IDs and dimensions; there is no fixture-name
whitelist and no suppression of numerical refusals. One new fixed-work charge
and the existing numeric decisions are charged to the caller's operation.

For P5 with `r=1/2`, the rail radius is `12.5 mm`. Every mouth's nearest point is
`17-4=13 mm`, so the finite cutters are clear despite their support intersections.
This edited public recipe now compiles and answers exact interior, exterior,
boundary and named-fillet queries. The original `r=2` has rail radius `14 mm`;
its mouths intrude, and this proof deliberately cannot admit it.

## Remaining first-step blocker

The authored P5 needs four rolling-envelope patches between the torus pieces,
where the ball rests on a mouth rim instead of on the flange plane. This is
real contact continuation, not a false-positive support intersection. The
existing B-rep implementation in `catallax-features::fillet_envelope` constructs
that topology, checked by `rim_contact_fillets`. Reusing its finished B-rep as
an allegedly independent native field would not close M18.

The field currently represents named boundaries with analytic support
polynomials and local trim/ray predicates. P5 continuation additionally needs a
native region representation with exhaustive coverage of the finite envelope
charts, their transitions, the mouth arcs and their material-side incidence.
`RollingEnvelope::side_at` alone is insufficient: its contract requires the
query's projection to lie on the finite chart, and it can refuse within the
chart enclosure width. A fitted carrier's position bound cannot be promoted to
an exact signed distance or an exact definition-side certificate. No such
promotion or blanket fallback is added here.

Consequently the authored P5, the tangent-mouth case `r=1`, and edits `r=5/2`
and `r=7/2` retain `fillet.unsupported-overflow`. Zero radius retains
`fillet.non-positive-radius`. The next required proof is the native finite
rolling-envelope solid region, including named transitions and unresolved-band
handling; this note is not a claim that the geometry is impossible.

Section 4.17 smooth blend is **not attempted after this unresolved first exit**.
It must read true operand distances; feeding existing Boolean residuals into
its circular formula would change the recipe's set. Wider field sweeps/lofts,
complete P1–P5 parity and production parallel compilation remain open.

## Executed validation

The retained [before probe](evidence/m18-clearance-before.txt) compiled the public
P5 recipe with the listed radius edits under 20M before the change; all five
refused overflow. That probe measured compilation only. The
[final matrix](evidence/m18-clearance-matrix.txt) includes compilation and six
point queries plus a named-face query for each successful row, so its work is
not a like-for-like speed comparison with the before probe.

| Case | Final outcome | Work per canonical job |
| --- | --- | ---: |
| P5 `r=1/2` | Builds; six analytic membership controls and named boundary pass | 11283 |
| Same, translated `(7,11,-13)` mm | Same expected geometry classifications and name membership | 11283 |
| Same, rotated 90° about z | Same expected geometry classifications and name membership | 10359 |
| P5 `r=1` | `fillet.unsupported-overflow` | 3298 |
| Authored P5 `r=2` | `fillet.unsupported-overflow` | 3367 |
| P5 `r=5/2` | `fillet.unsupported-overflow` | 3367 |
| P5 `r=7/2` | `fillet.unsupported-overflow` | 2811 |
| P5 `r=0` | `fillet.non-positive-radius` | 784 |

The control meridian uses direction `(3/5,4/5)` about the housing axis. Its
fillet centre is at radial/axial coordinates `(12.5,5.5)` mm, and the boundary
control `(12.2,5.1)` is exactly radius `0.5` away by the 3–4–5 identity. Nearby
inside/outside controls, the bore, a countersink interior and ordinary housing
material are checked separately. This is analytic control evidence, not sampled
global B-rep parity or an export check.

The validation workload uses eight fixed jobs. `OperationAllowance` reserves
8M once inside the unchanged 20M request and transfers 1M per job without
refunds or replacement allowances. Actual scoped workers are 1, 2, 4 and the
default available count (4 here), each nonempty. Records return in canonical
case order. Equality compares all `Query` fields, full `FaceQuery`, ordered
names, complete `Refusal` including attached undecided evidence, work and
exhaustion. A second run caps each job at 100 units; all eight refuse with equal
payloads/work across the same worker matrix. This is real parallel qualification
of the validation batch, not a claim that production recipe compilation now
uses all cores or that finite tests prove every thread count.

Executed commands:

```sh
cargo test -p catallax-field --test m18_clearance --test fillets
cargo clippy -p catallax-field --all-targets -- -D warnings
cargo run --release -p catallax-field --example m18_probe
cargo fmt --all --check
git diff --check
```

Seven tests pass: the new matrix control and six existing fillet regressions,
including P3, a convex ring, reflected/patterned fillets and multi-face corners.
Strict Clippy, formatting and diff checks pass. Cargo emits the existing
duplicate-bin manifest warning. No C++ oracle, protected corpus, saved inventory
edit, tolerance/allowance change, hosted CI launch, new production proof route
for an envelope, or smooth-blend acceptance is included. Hosted integration and
preservation checks remain unverified.
