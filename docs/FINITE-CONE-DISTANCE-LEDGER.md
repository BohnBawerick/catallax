# Finite conical face distance

Assignment: Firstmate plan-3 lane 5, correlation `0342c6dd51d90daf`.
Base: `1089120e14f777a082c9f09308c10dd81d90a74f`.
Predecessor #114 is included at merge
`a3f95cb400016416e8261f93ad7dbccf37804770` (source head
`82a256e09566fb7ae871d06141b9d151bdeeaa88`).

Only the shapes analytic distance implementation, measurement distance adapter,
focused tests/query benchmark, public synthetic fixtures and this lane's documents
are changed. No private inputs or held-out qualification are used.

## Defect and proof

Before this change, complete perpendicular-capped cones and frusta were not
recognized by `Finite`. `PreparedFace` used parameter subdivision even though
their complete nearest-point problem is a line-segment problem. Baseline and
candidate measurements use the identical public query harness, with the real
`PreparedFace` boundary-excluding fallback on baseline.

Let `a` be the unit axis, `k > 0` the slope, `z = (p-apex)·a` and `rho` the
perpendicular query distance. At any admitted height `h >= 0`, the surface has
radius `k*h`. Its azimuthal squared distance is
`rho² + (k*h)² - 2*rho*k*h*cos(theta) + (z-h)²`.
Since `rho*k*h >= 0`, aligned azimuth attains the minimum; at the axis or apex
all azimuths tie. Thus the complete minimum is the strictly convex quadratic
`(rho-k*h)² + (z-h)²` on the retained interval `[lo,hi]`. Its unique height is
`clamp((z+k*rho)/(1+k²),lo,hi)`. This includes the two rims, axis ties and an
explicitly retained apex. There is no parameter subdivision or angular sampling.
The lower bound encloses this exact minimum; the upper bound measures distance
to its exact admitted witness using the existing budgeted enclosure routine.

Recognition checks all boundary loops, positive-nappe circle incidence and exact
coaxial geometry. With outward normal `r-k*a` and natural circle tangent `a×r`,
`n×t` points along `a+k*r`, into increasing height. Combining the circle normal
with effective face/edge/coedge senses must yield one lower and one upper bound.
A vertex loop can supply only the exact apex lower bound. An excluded angular
sector cannot pass the complete-ring inventory check. The opposite algebraic
nappe fails positive-height recognition; witnesses are constrained to `[lo,hi]`,
so carrier continuations beyond either rim cannot donate points.

## Case inventory and retained limits

Tests cover interior projection, both rims and endpoint ties; axis ties, apex
and the opposite nappe; irrational azimuths; exact rigid poses and scales;
reversed cone axes, growing/shrinking stock frusta, all eight effective sense
combinations; invalid incidence/nappe/orientation and collapsed-loop inventories;
angular-sector fallback; invalid controls, zero work and insufficient-width
fallback equivalence. Independent exact rational and directed MPFR oracles
check complete result enclosures, not just overlapping intervals. Witness checks
use finite membership and a separate positive-nappe carrier equation.

General trims retain their complete fallback. Query support does not enable any
new shell construction. Zero query cells is not complete-click speed evidence.
All original defaults, tolerances, width, work, cell and vertex caps are retained.
The separate owner-arranged approved477 qualification is not run in this lane.

## Validation evidence

The first focused runtime run found a mistaken new-test expectation: a one-cell
fallback cannot initialize the complete multi-chart cone cover and returns zero
visited cells, not one. The full fallback result was already byte-for-byte equal
with/without the direct attempt. The regression now asserts the existing zero-cell
nonconverged behavior with the same one-cell cap; neither production behavior nor
a control was weakened. That original failed run is preserved with the evidence.
The subsequent expanded focused run passed 40 distance unit tests plus three filtered
integration tests. The existing spherical, cylinder and elliptic regressions pass.
Full-workspace and public artifact checks are reported with exact source identities in the draft PR; this ledger records the reproducible query experiment below.


## Matched public query experiment

Build `cargo bench -p catallax-measure --bench distance --no-run` on both source
revisions using the exact same benchmark file. Run the preserved executables
with `CATALLAX_CONICAL_DISTANCE_BENCH=1`. The baseline runtime is the named base
above; only its development harness is copied from the candidate. One warmup
and three alternating measured pairs ran as fresh processes pinned to one core.
Each query uses a fresh Decider. Per-pair load averages and raw rational bounds
are preserved outside the repository. Exact rational squared-distance oracles
check every complete enclosure. Counts and exact bounds repeat deterministically.

| Query | Baseline cells / query work | Baseline query ms | Direct query work | Direct query ms |
| --- | ---: | ---: | ---: | ---: |
| frustum-interior | 0 / 7,982 | 0.123 | 2,847 | 0.042 |
| frustum-lower-rim | 96 / 437,896 | 13.054 | 3,071 | 0.031 |
| frustum-upper-rim | 410 / 1,381,444 | 47.147 | 3,070 | 0.029 |
| frustum-axis | 0 / 10,498 | 0.124 | 3,743 | 0.038 |
| frustum-opposite-nappe | 84 / 397,444 | 11.625 | 3,070 | 0.023 |
| tip-interior | 0 / 7,982 | 0.097 | 2,847 | 0.031 |
| tip-lower-rim | 104 / 462,073 | 12.804 | 3,071 | 0.028 |
| tip-upper-rim | 0 / 8,998 | 0.103 | 2,847 | 0.025 |
| tip-axis | 0 / 10,498 | 0.110 | 3,687 | 0.035 |
| tip-opposite-nappe | 0 / 8,998 | 0.098 | 2,847 | 0.025 |

Every direct query converges with zero cells. The baseline also converges in all
ten rows: this is a direct-route gap, not a new unsupported-input refusal.
The identical controls are width 0.001 mm, max_cells 4,096, and the existing
library numeric defaults max_work 1,000,000,000 / max_precision_bits 16,384.
These library controls are separate from the unchanged 1M scoreboard, 20M app
and 100M verification-request regimes; no general application speed gate is
claimed from these shared-cloud measurements.

Preparation is outside query clocks and work: frustum preparation increases
from 3,274 to 3,499 work units and tip preparation from 3,371 to 3,590. Model
construction has its own clock and is excluded from both. Preparation is about
0.4–1.0 ms in this experiment. Fresh processes, preparation, construction and
query work are reported separately; zero cells alone is not whole-click proof.

Commands required for completion include release focused distance tests, the
full workspace release tests and doctests, workspace/all-target Clippy with
warnings denied, formatting, real verification-driver and independent STEP
regressions, edge-hunt/release-policy tests and the generated-page check.
Representative public strict comparisons retain manifests, original requests,
work/operation logs and actual STEP/OBJ/STL bytes. Fresh 1/2/4-worker and cold/warm
runs check determinism and validated artifact reuse separately. Any inherited
failure remains a failure and is reported in the PR; no allowance, assertion or
process cap is relaxed to pass it.

## Whole-ring property integration

The drill-point measurement path exposed a separate admission gap in the existing
closed-circle antiderivatives. A cone's apex coordinates can retain algebraic
construction expressions even when their transverse components are exactly zero.
The squared coaxiality expression then has no cached rational value. Treating that
cache miss as an unsupported ring sent volume, area, centroid and inertia through
Taylor panels.

For an already complete circular boundary with a parallel normal, the adapter now
asks the existing Decider to prove
`|(circle.center - cone.apex) cross cone.axis|² = 0` when the rational cache cannot
settle it. Only a proved zero admits the existing ring formulas. A proved nonzero
or unresolved identity retains quadrature; a refusal propagates. This does not
change finite-cone distance evaluation, trim membership, integration formulas,
requested accuracy, or any work allowance.

The public diagnostic case is a 30 × 30 × 10 mm block with one centred,
diameter-5 mm, depth-6 mm hole and a 118-degree drill tip. At base
`7a347113d34131a54dc59af2c6f073da2af9099f`, profiling confirmed a complete ring
whose rational coaxiality check returned no value while the exact predicate
proved coincidence. The profile spent 35.3 seconds in integration attempts and
2.0 seconds evaluating their nonrational literals on a contended host. A literal
cache therefore was not the selected fix.

The uninstrumented candidate's full mass-properties probe, including inertia,
completed in 599 ms. All 14 returned enclosure entries were contained in the
baseline enclosures. These are diagnostic observations, not SolidDraft click or
release timing acceptance. The application scenario remains the speed gate and
must use its ordinary 20-million-work policy.

The regression suite includes a complete cone translated by algebraic expressions
that equal zero, checked against closed-form mass properties and face centroids.
Separate positive and negative offsets of `2^-40` must not be admitted as coaxial.
These checks cover the shared admission condition for area, second moments and
face-centroid integrals without rounding a geometric predicate.
