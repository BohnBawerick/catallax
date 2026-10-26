# P1 Boolean proof costs, 26 October

This is a construction-cost investigation against `b1eb6f3c5cf454de008f8aef5162701d1dc77070`, not full-P1 acceptance. The complete original recipe and the configured boss-radius edit remain unchanged. Every construction request has the same 1,000,000-unit allowance.

## Baseline measurement

The [per-stage CSV](2026-10-26-p1-baseline-stages.csv) records all 82 reached Boolean stages across the original and the radius-9 edit, including both refused stages. Temporary development probes snapshot `OperationSpan` counters and elapsed time before propagating a refusal. They do not rely on `BooleanAnswer.report`, which contains only successful stages and is discarded by feature evaluation. No profiling probes remain in production code.

`work` is accepted request work, not the decider's cumulative attempted reservations. Percentages use the stage's decision count. Exact nonzero and exact zero decisions are separate columns. A zero-decision `Keep` stage has zero shares. Unreached stages have no row, rather than an invented zero cost.

| Baseline revision | Accepted construction work | Fact decisions | Filter decisions | Refined decisions | Exact nonzero decisions | Exact zero decisions |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Original | 999,999 | 4.871% | 75.205% | 6.231% | 4.588% | 9.103% |
| Boss radius 9 | 999,997 | 4.877% | 75.225% | 6.184% | 4.595% | 9.118% |

Each row also has one refused decision. There are no separation or not-exact decisions in these runs.

The original spends 279,884 work in `slot-cut-0` intersection and another 427,876 in the refused `slot-cut-1` intersection. The edit spends 275,831 and 431,927. Curve/surface intersections account for 230,673 and 404,027 work in those original stages. Whole-support intersections outside the finite boss and slot faces create expensive, irrelevant section curves.

The proposed planar endpoint-incidence reuse accounts for only 1,035 work across the reached chain, about 0.10% of the allowance. It is not selected as the first optimization. Its incidence fact would establish only carrier membership, not trim membership or orientation.

## Compatibility-preserving routes

- Use the circle/cylinder radical-plane calculation only to prove carrier separation or complete containment when the circle's normal is parallel to the cylinder axis. Offset axes are covered; coaxiality is not assumed. Tangency, overlap and unsupported orientations retain the original nonempty-hit solver and its RUR point expressions.
- Keep the original support-pair enumeration and rational-target line route. Omitting auxiliary support queries or replacing exact point expressions changed previously reported measurement enclosures, including wider bounds. Mathematical equivalence of the carriers does not guarantee identical interval arithmetic.
- Reuse the existing prepared signature workspace across one imprinted face's probes, with the same append-only model and decider. Cached projection proofs do not reuse a winding result for another point.

## Current construction smoke

Both unchanged recipes still refuse at `slot-cut-1`. The matched prefix before that Boolean is slightly cheaper; this is not completed construction or progress to a later feature.

| Recipe | Baseline work before `slot-cut-1` | Current work before `slot-cut-1` | Current accepted request work |
| --- | ---: | ---: | ---: |
| Original | 572,123 | 571,334 | 999,996 |
| Boss radius 9 | 568,070 | 567,281 | 999,996 |

A larger refused-operation count can mean that earlier work became cheaper and the operation reached more queries. It is not a regression measurement on equal completed work. The baseline CSV remains the per-stage profile; the prefix columns above come from the real CLI operation log, not a newly inferred Boolean-stage report.

## Verification and limits

The release intersection and Boolean suites pass 177 tests. One existing long-run qualification covering 18,392 Booleans remains ignored. The added regressions check exact roots, order and contacts for offset circle/cylinder sections, including internal/external tangency, nesting, containment and transformed poses; algebraic-line regressions check crossing, tangency and separation against rational cylinders.

All timings in the CSV are development observations on a shared host, pinned to one allowed core. They are not matched-pair latency qualification. Repeated observed construction waits vary substantially with ambient load and warmup. No speed acceptance follows from the timing columns.

P1 is not yet a completed construction or edit at 1M. Its unproved `preview_body` delivery and background proof remain unreached, so preview wait and proof time are unmeasured, not zero. The [26 October preview/proof contract](../ARCHITECTURE.md#preview-and-deferred-proof-decision-of-2026-10-26) remains unchanged. No allowance, tolerance, output accuracy or geometric definition is relaxed.

The real `catallax-verify` CLI replays both complete portable questions with their retained 1M allowances. Both honestly refuse, without producing a preview or an export. The exercised compatibility-preserving binary's SHA-256 is `d5e6355a2ce77cc115f00a27cc925dc82792a56332a5d8cfd5eadf98784ae531`.

## Public-case compatibility

The 41 reported public cases were replayed against the baseline and the compatibility-preserving variant with the same 100M verification-request allowance. This is a separate contract from the 1M P1 construction allowance.

Every prior answered manifest is identical after excluding kernel identity and statistics: topology, exact property endpoints, marks, STEP metadata and question identity are unchanged. One previous budget refusal now answers; all other refusals are unchanged. Both runs have zero crashes and zero cache hits. The earlier full point-solving and support-pruning experiments are not part of this variant.

All 35 resulting STEP files pass the calibrated standard-reader validity check with expected face and solid counts. The 34 pre-existing STEP files are byte-identical to baseline. This uses standard reader processing and is not a raw no-healing validity claim.

The geometry regressions and public replay exercise both the new exclusions and the unchanged nonempty-hit fallbacks. Full-corpus qualification remains separate from these 41 public cases.
