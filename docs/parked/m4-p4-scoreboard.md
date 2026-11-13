# P4 scoreboard work verdict

The remaining P4 work is parked with a measured `numbers.operation-budget`
refusal. The complete canonical request needs **1,072,166 work units**, 72,166
above the unchanged 1,000,000-unit scoreboard allowance. Exact construction
now fits at 917,132 units. At 1M the request completes construction and summary,
then refuses during STEP export. It leaves no STEP file. M4 remains open.

This records the remaining-work disposition for joint release review. It does
not change the M4 must-show list, certify a refused result or raise an allowance.
Resuming requires a later milestone to un-park the work under
[PLAN section 10.1](../PLAN.md#101-milestones).

## What main already delivered

The starting main revision was `015f64b8cf3b34626c28ab0a96cb9e2393ec51c3`.
It included the [cylindrical probe proofs](https://github.com/BohnBawerick/catallax/pull/481)
and [request enclosure reuse](https://github.com/BohnBawerick/catallax/pull/540),
which left the complete request at 1,366,808 units.

Main subsequently landed the [P1 proof-cost work](https://github.com/BohnBawerick/catallax/pull/546).
Its shared arithmetic and geometry changes also reduce P4 by 294,642 units.
This report retains that work and records its effect; it does not duplicate
those changes. Main also already contains the span-observer borrow repair and
its regression test, so this PR adds no kernel code. Recipes, tolerances,
allowances, comparators and the 477-entry saved-answer inventory are unchanged.

The measured source is `174f3474cbb3bca69a10694e9ed19137e4020f36`, rebased on
main `6fa7896a`. [Portable evidence](../evidence/m4-p4-remaining.json) records its
binary hash, fresh outcomes, caller inventory, stage partition and artifact
hashes. The [cloud validation record](https://github.com/BohnBawerick/catallax/tree/cloud/run-cx-m4-p4-remaining--p4-current-main)
records 95 passing numbers tests with `work-span`, both successful release
builds and the completed collector.

## Where the request stands

| Operation | Before the main refresh | Current source |
| --- | ---: | ---: |
| Shell | 62,617 | 60,290 |
| Rounding Boolean | 1,141,773 | 856,842 |
| Exact construction | 1,204,390 | 917,132 |
| Summary | 32 | 32 |
| STEP export | 162,386 | 155,002 |
| Complete request | 1,366,808 | 1,072,166 |

Construction includes shell and rounding Boolean. At 1M, the operation log
records both construction operations and summary as answered. The request
then spends its remaining 82,836 units in STEP and refuses at 1,000,000.
The manifest remains refused, with no returned result or STEP artifact. The
higher-cap answer has one solid, 15 faces, 36 edges and 24 vertices.

The independent span observer partitions the completed request as follows.
Its totals agree with the work observer and production at both 1M and 100M;
there is zero unattributed work.

| Completed-request stage | Work |
| --- | ---: |
| Intersect | 420,285 |
| Imprint | 409,564 |
| Classify | 83,205 |
| Sew | 244 |
| Other construction | 3,834 |
| Summary | 32 |
| STEP export | 155,002 |
| Total | 1,072,166 |

These stage rows include the shell's internal Boolean. The final refusal
location does not mean all required savings must come from export. Earlier
savings would also leave more of the shared request allowance for STEP.

## The residual inventory

The earlier enclosure-reuse investigation measured 1,043,962 baseline units
outside three leading proof rows. Enclosure reuse reduced that residual to
896,550. The current source reduces it again to **717,835**. Those older
baseline units are not today's remaining cost.

| Charged proof row | Work |
| --- | ---: |
| Probe coefficient separated enclosures | 89,693 |
| Circular-chart constant enclosures | 154,203 |
| Root isolation at the intersection solver | 110,435 |
| All other charged decisions | 717,835 |
| Complete request | 1,072,166 |

The three known rows total 354,331. Keeping all other work fixed would require
them to fit at most 282,165 units, a further 72,166-unit reduction. The old
argument that deleting those three rows could not close the gap no longer
applies after the landed reductions.

| Leading residual group | Work |
| --- | ---: |
| Stored point/curve incidence decisions, `store.rs:851` | 81,063 |
| Stored point equality decisions, `store.rs:948` | 102,317 |
| Fast face-region decisions, `fast_region.rs` | 123,616 |
| STEP number enclosures, `number.rs:130` | 55,571 |
| STEP parameter-curve enclosures, `pcurve.rs:730` | 49,671 |
| All remaining residual rows | 305,597 |
| Residual total | 717,835 |

These are exclusive numeric-caller charges, not whole-stage costs. Nested
checkpoint calls are charged separately. The known probe row includes only
its separated-enclosure label; other decisions at that source location belong
to the residual. The portable evidence retains every charged caller and label.

## What resuming must prove

Keeping construction and summary fixed leaves an export budget of **82,836**
units, down from 155,002. The two leading STEP enclosure rows together cost
105,242. If every other export charge stayed fixed, those rows would have to
fall to at most 33,076. Eliminating either row alone, even at zero cost, would
still leave the request above 1M. These are arithmetic bounds on the measured
partition, not implemented optimizations or lower bounds on future algorithms.

Keeping export and summary fixed instead requires construction at or below
844,966 units. Keeping shell fixed as well requires the rounding Boolean at
or below 784,676, down from 856,842. Shell work alone cannot close the gap:
even removing all 60,290 units would leave 1,011,876.

No further certified replacement was established here. The retained calls
prove geometric signs, incidence, whole-path clearance and the error bounds
on exported numbers. Replacing exact coincidence with tolerance, moving STEP
outside the counted request or removing cancellation charges would change the
contract rather than satisfy it. More threads do not reduce deterministic work.

Resume from main's retained proof reuse and inspect both export enclosure
callers as well as earlier construction. Measure the complete request after
every change. Keep the original canonical recipe, output bounds and STEP bytes.
The exit remains an exact completed answer at 1M, with the thickness-1.5 edit
reported separately and fresh agreement at 1, 2 and 4 workers. No higher-cap
answer counts as a 1M pass. This verdict does not assert that a future exact
algorithm cannot reach the allowance.

## Production checks

| Request | Work used | Outcome |
| --- | ---: | --- |
| P4 at 20M and 100M | 1,072,166 | Exact answer, original STEP hash |
| P4 thickness 1.5 at 100M | 1,091,184 | Exact answer |
| P4 at 1M | 1,000,000 | `numbers.operation-budget` in STEP, no STEP file |
| P4 thickness 1.5 at 1M | 1,000,000 | `numbers.operation-budget` in STEP, no STEP file |
| P4 at 100k | 99,998 | `numbers.operation-budget`, no STEP file |
| Near-miss control at 1M | 5,268 | Exact answer |

The 24-request matrix uses 1, 2 and 4 kernel threads with 1, 2 and 4 persistent
workers. Each configuration has four high-cap and four low-cap requests, so
every requested worker starts. Complete responses with clocks removed agree,
including work and operation logs. Strict STEP comparisons pass, all workers
exit zero and every request is fresh. Refusals leave no STEP file.

The P4 STEP SHA-256 remains
`f9a8b9a15a576a4989690f75e1bd9f697e052299b80b33c1bfb4fdc4634c7381`.
Answered manifests and STEP bytes match the preceding measured source
`e575b6ee7ff8c3a33657566cf0e5aad1d9837b7a`, excluding kernel commit identity
from manifests. Work counts and the refused operation change as reported above.
The protected 477-answer replay remains a separate CI gate. This focused matrix
does not claim the full release replay at 1/2/4 workers, field qualification or
SolidDraft acceptance. M26 timing is also separate. Raw evidence retained privately.

## Reproduction

The [collector](../../tools/research/m4_p4_verdict.py) runs the production
verifier through the existing driver, then the development work and span
observers separately. Set `PROOF` to an external evidence directory:

```sh
mkdir -p "$PROOF"
cargo build --release -p catallax-cli --bin catallax-verify -j 2
cp target/release/catallax-verify "$PROOF/catallax-verify"
cargo build --release -p catallax-cli --features work-span \
  --example work_profile --example work_span -j 2
python3 -B tools/research/m4_p4_verdict.py \
  --binary "$PROOF/catallax-verify" \
  --profile target/release/examples/work_profile \
  --span target/release/examples/work_span --out "$PROOF/run"
```

The collector checks actual STEP presence or absence, preserves per-request
allowances and requires the driver's strict comparisons to pass. It keeps
stderr separate from JSON and rejects an internal-error response. `portable.json`
contains no machine paths. Copying the production binary before the observer
build keeps production helper execution separate from the profiling feature.
