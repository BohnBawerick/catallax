# P4 scoreboard work verdict

The remaining P4 work is parked with a measured `numbers.operation-budget`
refusal. The complete canonical request needs 1,366,808 work units, 366,808
above the unchanged 1,000,000-unit scoreboard allowance. Exact construction
alone needs 1,204,390. P4 does not pass that gate, and M4 remains open.

This records the remaining-work disposition for the joint release review.
It does not change the M4 must-show list, certify a refused result, or raise
an allowance. Resuming this work requires a later milestone to un-park it
under [PLAN section 10.1](../PLAN.md#101-milestones).

## What already works on main

The starting main revision is
`015f64b8cf3b34626c28ab0a96cb9e2393ec51c3`.
It includes the [cylindrical probe proofs](https://github.com/BohnBawerick/catallax/pull/481)
and [request enclosure reuse](https://github.com/BohnBawerick/catallax/pull/540).
Those changes remain in place. Production arithmetic, recipes, tolerances,
work allowances, saved replay answers and comparators are unchanged.
The development-only span observer now releases the work counter's mutable
borrow before querying worker state. A real P4 observer run exposed the nested
borrow panic; a regression test checks admitted and refused charges. The
collector keeps stderr separate and rejects an internal-error response.

The earlier [enclosure-reuse profile](../M4-P4-ENCLOSURE-REUSE.md) measured
1,043,962 baseline units outside its three leading rows. Reuse removed
147,412 of those units, leaving 896,550. The fresh profile here starts after
that saving. The earlier baseline is not counted as today's remaining work.

## Remaining work

The observer accounts for every work unit, with zero unattributed work.
The following rows partition the complete request. The three known rows keep
the exact labels used in the earlier profile; other decisions at those source
locations belong to the residual rows.

| Charged proof | Work |
| --- | ---: |
| Probe coefficient separated enclosures | 205,620 |
| Circular-chart constant enclosures | 154,203 |
| Root isolation at the intersection solver | 110,435 |
| All other charged decisions | 896,550 |
| Complete request | 1,366,808 |

The independent span observer gives this stage partition, including the
shell's internal Boolean. Its total agrees with both the work observer and
the production verifier at 1M and 100M.

| Completed-request stage | Work |
| --- | ---: |
| Intersect | 462,478 |
| Imprint | 555,596 |
| Classify | 182,173 |
| Sew | 309 |
| Other construction | 3,834 |
| Summary | 32 |
| STEP export | 162,386 |
| Total | 1,366,808 |

Unlike the pre-reuse baseline, the residual now costs less than 1M. Reducing
only the three known rows could therefore close today's gap if their combined
cost fell from 470,258 to at most 103,450. Their old insufficiency argument
does not apply after enclosure reuse. Removing both coefficient and chart
enclosures entirely would still leave 1,006,985 units, so even that hypothetical
needs a further 6,985-unit saving in roots or residual work.

The residual work includes geometric signs, point/curve incidence, point
enclosures and STEP parameter curves. It is not an unmeasured overhead that
can be subtracted from the scoreboard.

| Leading residual group | Work |
| --- | ---: |
| Boolean geometric signs, nonzero | 170,695 |
| Boolean geometric signs, exact zero | 102,616 |
| Stored point/curve incidence signs | 143,448 |
| Stored point/curve incidence enclosure prefilter | 50,692 |
| STEP number enclosures at `number.rs:130` | 55,534 |
| STEP parameter-curve enclosures at `pcurve.rs:730` | 49,671 |
| All remaining residual rows | 323,894 |
| Residual total | 896,550 |

These are exclusive numeric-caller charges, not whole-stage costs. For
example, checkpoints inside a caller are attributed separately when nested
decider calls report them. The [portable evidence](../evidence/m4-p4-remaining.json)
retains the complete caller inventory and operation partition.

## What a reduction must accomplish

| Proposed scope, assuming its work vanished entirely | Maximum saving | Complete work left | Still over 1M |
| --- | ---: | ---: | ---: |
| Shell alone | 62,617 | 1,304,191 | 304,191 |
| Summary and STEP alone | 162,418 | 1,204,390 | 204,390 |
| Shell, summary and STEP together | 225,035 | 1,141,773 | 141,773 |
| Every exact-zero sign and every enclosure containing zero | 304,011 | 1,062,797 | 62,797 |

These are arithmetic ceilings for the recorded partition, not implemented
optimizations or lower bounds on future algorithms. They show why removing
only one of these cost classes cannot close the gap. In particular, replacing
exact coincidence with a tolerance decision would change the proof contract
and still would not remove enough recorded work, even at zero cost.

Keeping summary and export fixed requires construction at or below 837,582
units. Keeping shell fixed as well requires the rounding Boolean at or below
774,965, down from 1,141,773. A certified implementation therefore needs a
366,808-unit reduction across the Boolean proofs or a measured combination
of construction and export savings. Additional threads do not lower this
deterministic work count.

Source inspection found no remaining proof call that can simply be deleted:
geometric signs select sides and trim regions, incidence distinguishes curve
carriers, probe barriers certify whole-path clearance, and export enclosures
bound the written numbers. No further certified replacement was established
in this investigation. This is a bounded failure to finish the gate, not a
claim that a future exact algorithm cannot fit it.

## Production checks

The measured source is `e575b6ee7ff8c3a33657566cf0e5aad1d9837b7a`.
The [cloud validation record](https://github.com/BohnBawerick/catallax/tree/cloud/run-cx-m4-p4-remaining--p4-verdict-fixed)
records 91 passing numbers tests with `work-span`, both successful release
builds and the completed collector. Strict local Clippy passes with and without
the feature. The three focused span tests include the regression that first
failed with the original nested-borrow panic.

Fresh production requests reproduce the retained enclosure-reuse costs:

| Request | Work used | Outcome |
| --- | ---: | --- |
| P4 at 20M and 100M | 1,366,808 | Exact answer, original STEP hash |
| P4 thickness 1.5 at 100M | 1,444,166 | Exact answer |
| P4 at 1M | 999,993 | `numbers.operation-budget`, no STEP |
| P4 thickness 1.5 at 1M | 999,995 | `numbers.operation-budget`, no STEP |
| P4 at 100k | 99,984 | `numbers.operation-budget`, no STEP |
| Near-miss control at 1M | 6,980 | Exact answer |

The 24-request matrix uses 1, 2 and 4 kernel threads with 1, 2 and 4 persistent
workers. Each configuration has four high-cap and four low-cap requests, so
every requested worker starts. Complete responses with clocks removed agree,
including work and operation logs. Strict STEP comparisons pass, all workers
exit zero and every request is fresh. Refusals leave no STEP file.
All 24 production responses and STEP artifacts also match the preceding
source `9e1e01aff00a58b935c3331921d1285052af9a9f`, with only clocks and kernel
commit identity excluded from response comparison. The observer repair changes
neither production work nor geometry.

The P4 STEP SHA-256 remains
`f9a8b9a15a576a4989690f75e1bd9f697e052299b80b33c1bfb4fdc4634c7381`.
The protected 477-answer replay is a separate CI gate. This focused matrix
does not claim that full replay, field qualification or SolidDraft acceptance.

## Reproduction

The [collector](../../tools/research/m4_p4_verdict.py) runs the production
verifier through the existing driver, then runs the development work and span
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
allowances, compares complete production responses with only `wall_ms`
removed, and requires the driver's strict comparisons to pass. It requires
fresh artifacts and zero worker exits. `portable.json` contains no machine
paths. The production binary is copied before the observer build so helper
execution is tested independently of the profiling feature.

## Resuming and acceptance

Resume from the retained cylinder chart and enclosure reuse, without changing
the canonical P4 recipe. Measure every changed proof caller as well as complete
construction, summary and export. A reduction in one caller is insufficient
if it adds work elsewhere or shifts a refusal to a later operation.

The exit remains a complete exact answer at 1M, with the thickness-1.5 edit
reported separately, strict STEP bytes and the 477 saved replay answers
preserved, and fresh agreement at 1, 2 and 4 workers. Retain typed refusals and
the absence of STEP output below the required work. A higher-cap answer is
diagnostic evidence, not a 1M pass. SolidDraft interaction, field comparison
and M26 timing remain separate checks. Raw evidence retained privately.
