# P1 work profile and allowance verdict

Canonical P1 now fits the unchanged 20M application allowance. Its eager build
uses **15,264,436 work units** with 0.01 mm display accuracy, leaving 4,735,564
units spare. Exact construction uses 13,629,161 units. Verification through
summary and strict STEP uses 15,250,315 units. The normalized manifest and strict
STEP bytes match the original main baseline. Full M2 acceptance still needs the
separate SolidDraft lane.

The input is the unchanged [canonical P1 recipe](../parts/P1-bracket.json).
The existing work-profile observer attributes every charged unit. The
[portable evidence](evidence/m2-p1-work-profile.json) records source and binary
identities, operation costs, numeric callers and fresh comparisons. Raw evidence
retained privately. These are deterministic work counts, with no timing claim.

## Reductions and current verdict

The numeric and geometric reductions are measured in sequence:

| Source | Eager application work | Exact construction work |
| --- | ---: | ---: |
| Original main `4d51b5c3` | 41,647,298 | 39,863,121 |
| Parallel-cylinder generator proof `ccebbf65` | 33,754,179 | 31,970,002 |
| Reuse separating coordinate enclosures `142b3ab0` | 29,866,897 | 28,082,720 |
| Request enclosure proofs `5d6f6040` | 15,264,436 | 13,629,161 |

The last reduction removes 26,382,862 units from the original eager build,
or 63.3%. Its 20M and diagnostic 100M application runs have identical work and
counts: 19 faces, 42 edges, 20 vertices, 2,800 display positions and 5,620
triangles. The recipe, 0.01 mm display accuracy, 200,000-vertex display cap,
application allowance and numeric policies are unchanged.

The generator proof recognizes retained intersection branches of two distinct
parallel cylinders. It proves the common axis direction and nonzero radial
offset, then uses the stored intersection witness to select the generator.
The existing cylinder-chart polynomial proves whole-path clearance. Rotated,
translated, nonunit antiparallel axes and both algebraic witness branches have
regression coverage. Transverse and coincident cylinders remain outside this
proof.

The point-test reduction accepts an outward coordinate enclosure that excludes
zero as proof that points are apart. Repeating a sign proof for that same
coordinate adds no evidence. Unseparated coordinates still use the existing
exact-sign and distance-bound path.

Every recipe evaluation now opts into request-scoped enclosure reuse. A cache
entry contains a completed node's ball, its own cold arithmetic cost and its
working precision. A bounded query visits its expression DAG once, charges
each node's cold cost against the same per-decision limit, and reuses the
arithmetic already paid by the request. Cache keys separate working precision,
numeric policy and known-defined assumptions. A proof made under a domain
assumption cannot answer a checked-domain query.

The cache keeps its existing 8,192-entry bound, pins node identity, starts empty
for each request and reserved allowance, and is cold on helpers. Cache access
makes speculative records unrecordable, so the writer computes against its own
cache history. Regressions check shared-DAG charges at the exact cold limit and
one unit below it, policy separation, undefined square roots, and request scope.

The largest remaining eager-build callers are root isolation in `solve.rs`
(2,018,929 units), chart cancellation checkpoints (1,569,666), refined signs in
intersection math (1,179,762), curve-membership enclosures (1,028,021), and
whole-path coefficient enclosures (796,554). The final bridge cut now costs
4,344,198 units. These remaining charges fit inside the allowance.

The current application allowance is met. This result supplies no SolidDraft
interaction or timing qualification. The following sections retain the initial
profile and its failed verdict as historical evidence.

## Current production preservation

The fresh production matrix uses code
`79f6f7119c4319cd0a9f37d859aced359dea80e8`, built with the default release
profile and without either observer feature. Four requests run in each setting:
near-miss and line-tangent controls at 1M, plus canonical P1 at 100M and 20M.

| Scheduling | Persistent workers / kernel threads | P1 work at each allowance |
| --- | --- | ---: |
| Default | 1/1, 2/2, 4/4 | 15,250,315 |
| `CATALLAX_RUN_AHEAD=loops` | 1/1, 2/2, 4/4 | 15,250,315 |

All 24 fresh requests answer. Their complete responses agree across settings
after removing clocks and build identity, including work, numeric statistics
and operation logs. Strict STEP bytes agree, all workers exit zero, and there
are no driver cache hits or internal errors. Each result's normalized manifest
and STEP also match the original main answer. The formerly refused P1 20M
request matches the original completed 100M answer. The controls now cost 5,268
and 35,801 units.

The strict P1 STEP SHA-256 remains
`46d56e4250e9525fb026112b03519151322f08b5d9310bf99672557b61995f08`.

The `p1_allowance` application regression builds the unchanged recipe with
display under `DEFAULT_BUILD_WORK` at 1, 2 and 4 kernel threads. All three runs
use 15,264,436 units and compare equal in names and complete mesh data, including
per-face accuracy reports and boundary polylines. Its assertion also fixes the
default allowance at 20M; it does not supply a diagnostic override.

The direct five-stage Boolean test also passes at 1, 2 and 4 threads. All 16
intermediate results match normal construction in topology, face/edge/vertex
names and strict STEP bytes. Both paths use the same request enclosure
preference. Their per-operation counts are identical across thread counts;
`slot-cut-4` uses 19 fewer units in direct execution, and the other 15 Boolean
counts match normal execution. These are fresh per-prefix diagnostic budgets,
not the complete-request allowance proof supplied by `p1_allowance`.

The repaired span observer matches production work exactly at 20M, with every
reported critical path consistent with its span:

| Current stage | Work |
| --- | ---: |
| Intersect | 8,229,388 |
| Imprint | 4,198,496 |
| Classify | 1,192,373 |
| Sew | 8,194 |
| Other construction | 710 |
| Summary and STEP | 1,621,154 |
| Total | 15,250,315 |

Its normalized manifest and strict STEP bytes match both production and the
original main result. Enclosure reuse exposed recursive destruction of long
observer paths after the answer artifacts had been written. A debugger trace
identified repeated `Arc<Segment>` destruction; iterative destruction fixes
that development-only failure. A 100,000-segment regression runs on a 128 KiB
stack and retains a shared prefix. All 94 numbers tests pass with observer
features, including four span tests. Focused production and observer Clippy
pass with warnings denied.

The [hosted preservation run](https://github.com/BohnBawerick/catallax/actions/runs/43302466267)
replays all 477 protected records freshly, with no cache hits, failures or
verified corrections. All 477 answers match. The four partitions cover 86
intersection, 217 mixed, 89 subtraction and 85 union cases, with no duplicate
or missing indices. This tests GitHub's merge
`d88f548223136840a83a0f38159ba942b9be4b3d`, whose parents include the qualified
head `79f6f7119c4319cd0a9f37d859aced359dea80e8`. The portable evidence retains
that distinction, the hosted binary identity and the corpus/policy digests.
No saved replay record changes.

## What main already supplied

Main already has the verification `work_profile` observer and the eager
application `app_work_profile` observer. The latter landed with the rational
circle loft work in `908a9e50`. The full P1 name-edit gate is also separate
existing coverage, recorded in [the M2 note](milestones/M2.md).

The 44,824,339-unit P1 row in [RUN-AHEAD](RUN-AHEAD.md#7-speed-on-p0-to-p5)
measures an older complete verification request. It is neither an application
build nor the initial profile below. On main
`4d51b5c36bdad29632a448585290cfaa627098ba`, the entries had these costs:

| Entry | Allowance | Work used | Outcome |
| --- | ---: | ---: | --- |
| Application, exact body and topology check | 100,000,000 diagnostic | 39,863,121 | Built, 19 faces |
| Application, eager build and display | 100,000,000 diagnostic | 41,647,298 | Built, 19 faces, 42 edges, 20 vertices |
| Application, eager build and display | 20,000,000 | 19,999,999 | `application.work-limit`, `slot-cut-6` |
| Verification, through summary and STEP | 100,000,000 diagnostic | 41,501,362 | Answered, exact |
| Verification, through summary and STEP | 20,000,000 | 19,999,999 | `numbers.operation-budget`, `slot-cut-6` |

The diagnostic allowance measures a completed request. It does not change the
application default or turn its refusal into a pass. Exact construction alone
exceeds 20M by 19,863,121 units. Application display and naming add 1,784,177.
Verification summary and STEP add 1,638,241 instead. The refused verification
request leaves no STEP.

## Initial work profile

The application operation log gives the following complete construction costs.
Stock cylinders, blocks and extrusions cost no units in this recipe. The two
sketches together cost 710.

| Boolean | Work |
| --- | ---: |
| `base` | 39,263 |
| `bossed` | 19,278 |
| `bored` | 23,976 |
| `bolt-cut-0` | 25,991 |
| `bolt-cut-1` | 28,287 |
| `bolt-cut-2` | 30,663 |
| `bolt-cut-3` | 33,119 |
| `slot-cut-0` | 35,682 |
| `slot-cut-1` | 1,317,836 |
| `slot-cut-2` | 2,414,985 |
| `slot-cut-3` | 3,705,938 |
| `slot-cut-4` | 4,482,792 |
| `slot-cut-5` | 4,625,866 |
| `slot-cut-6` | 5,138,565 |
| `slot-cut-7` | 4,861,556 |
| `bracket` | 13,078,614 |
| All Booleans | 39,862,411 |

The overlapping slot cylinders and final bridge cut dominate the request. The
final `bracket` cut alone uses 31.4% of the eager build's work. The complete
construction finishes before display begins at the diagnostic allowance.

The repaired span observer independently partitions the complete verification
request. Its work matches the production count exactly:

| Stage | Work |
| --- | ---: |
| Intersect | 17,432,253 |
| Imprint | 16,397,625 |
| Classify | 6,022,948 |
| Sew | 9,585 |
| Other construction | 710 |
| Summary and STEP | 1,638,241 |
| Total | 41,501,362 |

The largest exclusive numeric callers in the eager build are:

| Caller and decision | Work | Share of eager build |
| --- | ---: | ---: |
| `probe_path.rs:36`, separated coefficient enclosures | 10,747,414 | 25.8% |
| `store.rs:1157`, point-test enclosures | 4,028,115 | 9.7% |
| `store.rs:1172`, refined point-test signs | 3,887,282 | 9.3% |
| `geometry.rs:43`, refined signs | 2,967,989 | 7.1% |
| `solve.rs:268`, root isolation | 2,018,929 | 4.8% |
| `math.rs:65`, refined signs | 1,889,060 | 4.5% |

Caller locations refer to the measured source. Enclosures that separate zero
and nonzero signs account for most work. Exact-zero signs cost 2,640,458 units.
The coefficient caller already stops refinement when an enclosure proves its
sign. Its whole-path clearance proof cannot be replaced by sampled clearance.

## Initial failed verdict

The eager build needs a reduction of **21,647,298 units**, or 52.0%, to meet
20M. Removing all display and naming work would still leave 39,863,121 units.
Even a zero-cost lower bound for all three leading numeric callers leaves
22,984,487 units. A reduction confined to those callers cannot close the gate
with every other charge fixed.

The initial profile established no construction shortcut or numeric proof
reduction. Further work must reduce proofs across the overlapping slot cuts and
the final bridge cut while retaining their incidence, branch and whole-path
witnesses. Combining tools or replacing the slot definition must also preserve
every intermediate answer and its names. More threads do not reduce the charge.
The present result is a measured `application.work-limit` refusal, not a claim
that the 20M gate is mathematically impossible.

## Construction routes and reproduction

The opt-in [P1 qualification test](../crates/catallax-features/tests/p1_boolean_routes.rs)
compares every Boolean with direct execution of intersect, imprint, classify,
keep and sew on the same operand prefix. This bypasses the Boolean entry's
separated-union and apart-body construction shortcuts. Exact numeric predicates
inside the five stages remain active. The test checks topology, all face, edge
and vertex names, and strict STEP bytes for all 16 intermediate results.
Its 100M per-operation diagnostic allowance does not qualify a 20M build.

```sh
cargo test --release -p catallax-features --test p1_boolean_routes \
  -- --ignored --nocapture
cargo test --release -p catallax-numbers --features work-span --lib \
  work_span::tests -- --test-threads=1
cargo build --release -p catallax-cli --bin catallax-verify \
  --features work-span --example work_profile --example work_span \
  --example app_work_profile -j 2
```

Use an external directory for `PROOF`. These calls retain the original recipe
and distinguish the two request entries:

```sh
INPUT=parts/P1-bracket.json
CATALLAX_THREADS=1 target/release/examples/app_work_profile "$INPUT" \
  --max-work 20000000
CATALLAX_THREADS=1 target/release/examples/app_work_profile "$INPUT" \
  --max-work 100000000
CATALLAX_THREADS=1 target/release/examples/app_work_profile "$INPUT" \
  --max-work 100000000 --exact-only
CATALLAX_THREADS=1 target/release/examples/work_profile "$INPUT" \
  --out "$PROOF/profile" --max-work 100000000
CATALLAX_THREADS=1 target/release/examples/work_span "$INPUT" \
  --out "$PROOF/span" --max-work 100000000
```

The span observer initially returned `cli.internal-error` after one unit because
its budget-charge callback reread a mutably borrowed request budget. The callback
now runs after that borrow ends. A regression test checks accepted work and a
reservation that does not fit. The fix is confined to development observation;
production arithmetic, tolerances, allowances and saved replay records are
unchanged. All three span unit tests pass. A complete P1 observation accounts
for every production work unit and preserves the normalized manifest and strict
STEP bytes. The eager and exact-only application profiles also reproduce their
baseline outcomes and counts. These development observers use release
optimization level 1; the separate production binaries use the default release
profile. Neither comparison is a timing qualification.

## Initial observer preservation

The baseline is the main snapshot above. Candidate code
`3f07d79e48e71b916b570833598eedc285ccbcdf` changes development observation only.
Both production verifiers were built without either observer feature. Profiling
builds disable helper execution and cannot prove production thread behavior.

The fresh matrix has four ordered requests: P1 at 100M and 20M, plus the unchanged
near-miss and line-tangent controls at 1M. The controls let all four persistent
workers start. The candidate runs each request in these configurations:

| Boolean scheduling | Persistent workers / kernel threads | Result |
| --- | --- | --- |
| Default | 1/1, 2/2, 4/4 | Equal to baseline |
| `CATALLAX_RUN_AHEAD=loops` | 1/1, 2/2, 4/4 | Equal to baseline |

All 24 candidate requests match the baseline's complete responses after removing
clocks and build identity. Work, statistics, operation logs, normalized manifests
and strict STEP bytes agree. Every worker exits zero, with no cache hits or
internal errors. Every refused request has no STEP file. The baseline work-profile
manifest and STEP also match the production binary exactly.

The construction-route test passes at 1, 2 and 4 kernel threads. All 16 Boolean
results match direct five-stage execution in topology, names and strict STEP
bytes. Each normal operation also charges the same work as direct execution,
and those per-operation counts agree across threads. No P1 construction shortcut
is observed in this snapshot. These prefix checks start fresh diagnostic budgets;
their counts do not replace the complete-request profile above.

Focused production and span-feature Clippy pass with warnings denied.

The initial observer-only run did not cover SolidDraft interaction, independent
STEP import, field comparison or the full 477-case replay. It did not complete
M2. Current preservation results are recorded above.
