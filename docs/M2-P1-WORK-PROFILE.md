# P1 work profile and allowance verdict

Canonical P1 still exceeds the unchanged 20M application allowance. On main
`4d51b5c36bdad29632a448585290cfaa627098ba`, its eager application build uses
**41,647,298 work units** with 0.01 mm display accuracy. At 20M it refuses with
`application.work-limit` during `slot-cut-6`, after 19,999,999 accepted units.
M2 remains open.

The input is the unchanged [canonical P1 recipe](../parts/P1-bracket.json).
The existing work-profile observer attributes every charged unit. The
[portable evidence](evidence/m2-p1-work-profile.json) records source and binary
identities, operation costs, numeric callers and fresh comparisons. Raw evidence
retained privately. These are deterministic work counts, with no timing claim.

## What main already supplies

Main already has the verification `work_profile` observer and the eager
application `app_work_profile` observer. The latter landed with the rational
circle loft work in `908a9e50`. The full P1 name-edit gate is also separate
existing coverage, recorded in [the M2 note](milestones/M2.md).

The 44,824,339-unit P1 row in [RUN-AHEAD](RUN-AHEAD.md#7-speed-on-p0-to-p5)
measures an older complete verification request. It is neither the current
count nor an application build. The current entries have these costs:

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

## Where the work goes

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
| `store.rs:1157`, point-membership enclosures | 4,028,115 | 9.7% |
| `store.rs:1172`, refined point-membership signs | 3,887,282 | 9.3% |
| `geometry.rs:43`, refined signs | 2,967,989 | 7.1% |
| `solve.rs:268`, root isolation | 2,018,929 | 4.8% |
| `math.rs:65`, refined signs | 1,889,060 | 4.5% |

Caller locations refer to the measured source. Enclosures that separate zero
and nonzero signs account for most work. Exact-zero signs cost 2,640,458 units.
The coefficient caller already stops refinement when an enclosure proves its
sign. Its whole-path clearance proof cannot be replaced by sampled clearance.

## Verdict

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

## Production preservation

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

SolidDraft interaction, independent STEP import, field comparison and the full
477-case replay remain separate gates. This landing does not claim their result
or complete M2.
