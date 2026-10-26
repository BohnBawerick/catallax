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

## Selected routes

- Apply the existing conservative finite-face exclusion proofs before rational support-pair construction. Strictly separated boxes or a support enclosure excluding zero prove only this finite pair irrelevant. An overlapping box proves nothing; unsupported bounds retain the analytic route. The algebraic support-pair route is unchanged.
- Use the existing one-parameter line solver when the line itself has algebraic coordinates, even if the intersected surface is rational. Rational-line/rational-surface cases retain their existing solver. This preserves correlated coordinates instead of independently eliminating them.
- Intersect a circle perpendicular to a cylinder's axis through its exact radical plane and common chord. Parallel, offset axes are covered; coaxiality is not assumed. Signed squared chord height distinguishes two crossings, both tangent cases and disjoint carriers. Other orientations retain the general solver. Candidate branch, domain, contact and ordering checks still run.
- Reuse the existing prepared signature workspace across one imprinted face's probes, with the same append-only model and decider. Cached projection proofs do not reuse a winding result for another point.

## Measured construction progress

| Development snapshot | Original `slot-cut-0` work | Original `slot-cut-1` intersection work | First refusal |
| --- | ---: | ---: | --- |
| Baseline | 307,350 | 427,876 before refusal | `slot-cut-1`, Intersect |
| Finite-face pruning | 107,761 | 669,139 before refusal | `slot-cut-1`, Intersect |
| Plus algebraic-line routing | 107,761 | 392,395, completed | `slot-cut-1`, Imprint |
| Plus offset circular sections | 107,761 | 219,333, completed | `slot-cut-1`, Imprint |

A larger refused-stage count can mean that earlier work became cheaper and the stage reached more queries. It is not a regression measurement on equal completed work.

With prepared signature reuse, both complete recipes still refuse during `slot-cut-1` imprinting at 999,995 accepted work. In that run, signatures consume 313,107 work across the reached chain. Trace preparation is only 1,614 work. Almost all remaining signature cost is certified projection/subdivision, not trace construction. The two expensive first signatures consume 136,035 and 158,893 work. Reusing traces alone does not remove that proof cost.

## Verification and limits

The release intersection and Boolean suites pass 177 tests. One existing long-run qualification covering 18,392 Booleans remains ignored. The added regressions check exact roots, order and contacts for offset circle/cylinder sections, including internal/external tangency, nesting, containment and transformed poses; algebraic-line regressions check crossing, tangency and separation against rational cylinders.

All timings in the CSV are development observations on a shared host, pinned to one allowed core. They are not matched-pair latency qualification. Repeated observed construction waits vary substantially with ambient load and warmup. No speed acceptance follows from the timing columns.

P1 is not yet a completed construction or edit at 1M. Its unproved `preview_body` delivery and background proof remain unreached, so preview wait and proof time are unmeasured, not zero. The [26 October preview/proof contract](../ARCHITECTURE.md#preview-and-deferred-proof-decision-of-2026-10-26) remains unchanged. No allowance, tolerance, output accuracy or geometric definition is relaxed.

The real `catallax-verify` CLI, driven through the verification skill with both full inputs, also refuses at 999,995 work. The run has two cases, zero crashes and zero cache hits. Observed request waits are 742.425 ms and 734.291 ms; these are refusal waits, not preview latency. The exercised binary's SHA-256 is `82e0585c2d5e8b3631d688787eeb74de75c2d882f2835f0e21a184206064f0a1`.
