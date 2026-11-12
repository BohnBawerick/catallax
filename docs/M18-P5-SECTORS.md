# Authored P5: angular partition prerequisite, admission still blocked

Baseline: `cdc0a3884c59b2bea538c230e37639bd74cdb9a2` (11 November 2026).
The input is the public `parts/P5-flanged-housing.json`, radius **2**.
This increment does **not** admit that recipe. It adds an independently exact
angular-routing prerequisite to the merged finite rolling region, without
changing compiler acceptance or the face registry.

## What is proved

`rolling_region::Sectors::new` takes an explicit list of validated `Region`s.
It requires identical origin, oriented axis, body radius and ball radius.
For each pair it checks

`toward_i dot toward_j < cos(theta_i) cos(theta_j) - sin(theta_i) sin(theta_j)`.

Each half-angle is in `(0, pi/2)`, established by `Region::new`. Their sum
is in `(0, pi)`; cosine strictly decreases on this interval. The inequality
therefore proves strict separation of the two entire closed angular sectors,
including the periodic seam. Equality refuses; there is no tolerance. This
also means at most one sector can own any nonzero radial direction.

`locate` partitions all nonzero directions in the common plane into strict
mouth interiors, their exact start/end rays, and the complement labelled
`Torus`. It checks exact sector halfspaces and preserves unresolved numeric
comparisons as refusals. The label does not certify a material point or a
boundary. Invalid directions, incompatible rings and touching/overlapping
sectors refuse as `field.rolling-sectors-unresolved`.

`transition` returns the existing definition-native station at `u=-1` or `1`.
There, planar centre-to-mouth distance is exactly `m`, so centre height is
exactly `r`. The ordinary torus and finite mouth have the same rolling ball,
cylinder contact and plane/mouth contact at that endpoint. Tests check these
identities, not approximate coordinates. Endpoint agreement alone does **not**
prove agreement of the lower envelopes over their entire joining section.

All data are immutable after construction. Region indices retain input order;
they are not persistent recipe names. Construction charges every region and
pair, lookup every visited sector, and exact decisions use the caller's
decider and operation budget. Cancellation and numeric refusals propagate.
These short routing loops are serial. The existing expensive region interval
queries still use production `classify_batch`, all available cores by default,
with fixed per-input shares of one allowance.

## Exact remaining obstacle

P5's four cardinal mouths can pass this explicit angular partition, but that
list is not a compiler-derived incidence certificate. Before admission:

1. Discover every active mouth and adjacent material face from the recipe's
   Boolean ownership, not a hardcoded P5 signature. Account for all other
   source faces and finite tool extents; an angularly disjoint mouth does not
   exclude foreign-feature interaction with its rolling balls.
2. Prove the regularised material union specified by RECIPE 4.5.1, including
   the radius homotopy. The fixed-radius clipped lower-envelope set is not
   already a proof of that union. Prove the torus/finite-region joins over
   whole sections and rail/mouth arcs; coincident endpoint contacts are
   necessary but insufficient. The complement of angular sectors is only
   angular coverage, not material coverage.
3. Supply exact boundary queries for the envelope's minimizing station,
   including endpoints, multiple minima and unresolved cases. `Face` currently
   stores a supporting `Surface`, polynomial, gradient and trimmed domain;
   support and ray queries use those analytic objects. A fitted B-spline
   carrier, or the ordinary full torus, cannot stand in for this envelope.
   The full-torus counterexample at `(12.5,5.5,0)` in
   [the region note](M18-ROLLING-REGION.md) remains decisive.
4. Integrate those boundaries with ownership, persistent origins, selection,
   normals, rays and distance coverage. No unregistered envelope leaf may be
   admitted and then expose a named-subset distance as a whole-boundary answer.

None of these obligations is discharged by sampled membership, a fitted
carrier tolerance, or the new angular partition. The compiler therefore
retains `fillet.unsupported-overflow` for authored radius 2. The radius-1/2
edit is still admitted; the normal-cone exit and other existing refusals remain.

## Interface to M17

No new `Plan`, mass integration domain, face, name or distance authority is
provided. M17 must continue to report the P5 compilation refusal. `Sectors`
may later route candidate native boundary queries; its `Torus` and `Mouth`
labels must not be used to integrate volume or to infer a complete solid.
Future admission must deliver the ordinary checked `Plan` contract, including
complete boundary coverage, before M17 can use it. No mass algorithm changes
are included here.

## Focused reproduction

```text
cargo test --release -p catallax-field --test rolling_sectors --test rolling_region --test m18_clearance -- --nocapture
cargo clippy -p catallax-field --lib --test rolling_sectors --test rolling_region --test m18_clearance -- -D warnings
cargo fmt --all --check
git diff --check
```

The sector controls use four explicit P5-dimension mouths, translated controls,
both exact endpoint contacts, directions immediately to each side, complement
quadrants, overlapping/touching sectors, incompatible rings, invalid directions,
deadline and exhausted allowance. Existing region controls cover the finite
material/full-torus distinction, boundary refusal and normal-cone exit.
Existing P5 controls exercise actual compilation and named queries, including
the radius-1/2 edit under translation and rotation. Batch tests compare full
semantic records (including complete refusals), work and exhaustion at
1/2/3/4/default/9 workers. Nine is capped by the input count. These finite
controls do not prove arbitrary worker-count behavior; fixed index-based
allowance partitioning and ordered collection are the source invariants.

No timing, independent B-rep comparison, STEP export, hosted gate, saved-answer
replay, complete M18 exit or authored P5 admission is claimed.

## Executed evidence (11 November 2026)

The three focused release targets passed: `rolling_sectors` 3 tests,
`rolling_region` 6, `m18_clearance` 1. After adding the P5 work-count output,
`m18_clearance` and `rolling_sectors` were run again (4 passed).
Targeted Clippy with warnings denied, workspace formatting check and whitespace
check passed. The existing duplicate CLI binary-target manifest warning remains;
no CLI build or broad suite was requested.

Default workers were 4. At 1/2/3/4/default/9, the P5 matrix's per-job work was
`[11283,11283,10359,3298,3367,3367,2811,784]`, without exhaustion. Case order is
radius-1/2 original/translated/rotated, then radii `[1, 2, 5/2, 7/2, 0]` on the
original placement.
The eight fixed 1M child reservations charged 8M once against the unchanged
20M caller limit. At 100 per child, each exhausted at 100 and refused, with
800 charged total. The requested 9 workers started 8 nonempty workers.

The production region batch's per-input work at a 700,000 reservation was
`[184,156,156,6,156,2,2]` at every tested worker count, with one clipping-wall
refusal and no exhaustion. The caller charged 700,002 (reservation plus two
checks). At reservation 7 all seven jobs exhausted at 1; at 0 all exhausted
at 0. Requested 9 workers started 7. Tests compare full records, including
refusal evidence, rather than the abbreviated diagnostic rendering.

Final release test executable SHA-256 identities:

| Target | SHA-256 |
| --- | --- |
| `rolling_sectors` | `0e542da41a120c8b7a1cebb74ef5a2801fd0aab48f2092c83b354fbb203fbfe3` |
| `rolling_region` | `6ed114ff0f9d2d32f70284916b750427f8901d5ad0c23c9975c16976def92e4b` |
| `m18_clearance` | `c305b4e7776ad8003ed822db919a149387aedabd086b716d30c5df2710aded1a` |

| Requested acceptance | Result |
| --- | --- |
| Authored radius-2 P5 field admission | Still refused; not complete |
| Explicit finite-sector angular partition and transition endpoints | Implemented, exact-source argument and focused controls |
| Exhaustive recipe incidence, material homotopy and whole-section joins | Open |
| Native envelope names/rays/distance and M17 P5 domain | Open; no substitute authority exposed |
| Existing admitted edit, refusal and normal-cone controls | Preserved in focused execution |
| Worker/accounting and exhausted allowance comparisons | Equal at tested counts; no universal runtime claim |
