# R3 generated recipes with built-in checks

This Rust-only companion to the edge-hunt outcome census implements a bounded
first tranche of PLAN.md sections 12.3 and 12.6 R3. It creates native public-format
recipes and checks successful measurement enclosures against independent exact
mathematics. It neither reads a saved-answer inventory nor uses an external CAD
oracle. These finite checks are evidence, **not a proof that the kernel never
returns a wrong proven answer**. No production algorithms or allowances change.

## Standing-lane invocation

From the repository root:

```sh
cargo test --release -p catallax-features --example generated_recipe_checks -j 2
cargo build --release -p catallax-features --example generated_recipe_checks -j 2
# SEED START COUNT NEW_OUTPUT_DIRECTORY; choose a fresh scratch output directory.
timeout 600s target/release/examples/generated_recipe_checks 20261111 0 96 "$R3_OUTPUT"
```

The caller supplies `R3_OUTPUT`; the directory must not exist. A standing/nightly
runner can invoke this command and retain `results.json` counts by outcome. No
scheduler or duplicate CI workflow is installed by this tranche. A second fixed
seed is `1414213562`. Increase count explicitly, up to 10,000; counts never renew
inside an invocation. The maximum is a safety ceiling, not a recommended nightly
size. Each generated case runs two variants (original and translated), or three
for union/intersection (also reversed operands), at each of 1, 2 and 4 workers.
Workers run independent single-threaded requests, not parallel kernel arithmetic.
Each build plus measurement has a fresh **20,000,000-unit** allowance, matching the
ordinary application build default; default measurement accuracy is unchanged.
At most four requests execute concurrently. Use an outer wall timeout and the
host's ordinary process/memory limits; the work limit is not a wall-time or memory
guarantee. A killed run, missing report or nonzero exit is incomplete/failed,
never passing evidence. All recipes are written before any kernel execution.

The generator is versioned by source commit. SplitMix64 streams derive from the
seed and case index; larger runs preserve the prefix. Record the exact Git head,
command and exit status beside the output. To replay case N, use the same seed,
START=N and COUNT=1 with a new directory. Every `case-N-V.json` is also a complete
`catallax-recipe/3` input usable by the existing recipe verification interface.
Output includes rational recipe coordinates, parameters, classification,
measurement enclosures/marks or refusal detail, and deterministic arithmetic work.
It deliberately omits timing from worker comparisons.

## Independent checks and coverage

Four operations cycle: primitive block, union, intersection, subtraction. Six
placement strata cycle every 24 cases: identical operands, positive gap, face
contact, corner containment and two general placements. Sizes and placements
vary with the seed. Denominators cycle through 1, 2, 8 every 24 cases. Negative
positions, empty subtraction/intersection, disconnected unions and coplanar faces
are included. No hidden fixture or customer input is needed.

The reference partitions space at the two boxes' coordinate endpoints. At most
27 open rectangular cells exist. Each cell's membership comes from integer
coordinate comparisons and the Boolean truth table. With cell side lengths
`dx, dy, dz` and scale denominator `d`:

- volume is `sum(dx*dy*dz) / d^3`;
- area is the sum of exposed cell-face areas divided by `d^2`; internal faces
  contribute zero, including coplanar faces between adjacent cells;
- twice first moment on axis i is
  `sum(dx*dy*dz*(low_i+high_i)) / d^4`;
- centroid is first moment / volume, absent exactly when volume is zero.

Reference membership, integer sums and face adjacency do not use kernel geometry,
topology, integration or enclosure code. Only conversion of the final rational
ratio uses the shared numbers type; a defect in that arithmetic remains a shared
trust dependency. Unit controls independently assert a 2x3x4 block's volume 24,
area 52 and centroid (1, 3/2, 2), equal subtraction's empty result, face-joined
area 80 and separated area 104. Fault injections move volume/area outside their
expected enclosures and remove the centroid; each must fail its named check.

Each returned volume, area and centroid enclosure must contain the exact answer,
including marked answers. Translation by 7/d on every axis must preserve volume
and area and translate centroid; union/intersection operand reversal must preserve
these properties. These are property-level metamorphic checks against the same
independent reference, not an assertion of identical topology or naming across
equivalent recipes. At 1/2/4 workers the *same input's* full measurement debug
record (including marks), refusal detail, outcome and work must match exactly.

## Outcomes and reduction

`checked` means all requested analytic checks passed. Build/measurement typed
refusals are counted separately by stage and code. Invalid recipe, invalid
topology, panic, wrong proven volume/area/centroid, budget contract, worker mismatch
and coverage failure fail the invocation. Every operation present in the selected
range must have at least one successful base variant, so an all-refusal run cannot
pass. The 24-case smoke additionally requires every base variant to succeed.

For each analytic/structural failing variant, up to twelve greedy reductions try
zeroing placement components and reducing side lengths to one. A reduction is
kept only if it reproduces the same failure class under the same allowance. The
final `reduced-N-V.json` and its parameter row are retained beside the original.
This is a bounded smaller reproducer, **not a globally minimal recipe**. Worker
mismatches retain both observations and the original inputs; they are not shrunk.

Limits: axis-aligned boxes only, no curved surfaces or relative rotation; no
point-membership proof, export check or general independent body checker. Correct
volume, area and centroid do not uniquely determine a solid. No saved-answer replay,
C++ oracle, target-part qualification or full hosted gate is implied. Future
tranches can add independently justified analytic families and deeper reductions.
