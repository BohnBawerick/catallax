# Hard intersections

Protocol recorded before prototype code. Work bound: one measured round of both candidates against the baseline, with at most one repeat and only if the measurement itself was wrong. The verdict lands before the hard-intersection repair work it informs starts. This is an experiment, not a change to the kernel's representation.

## Question and candidates

Can either representation avoid the large-number polynomial refusals in equal-radius and seam intersections without changing an earlier proven answer?

- **Implicit:** retain the defining surface pair and component selector; decide membership, ordering and separation against those definitions without expanding a global elimination polynomial wherever a direct certificate suffices.
- **Certified approximation:** retain the same definition beside certified local curve enclosures; use those enclosures only when they prove the needed decision. Refine or use exact geometry when they cannot. Target deviation is 100 nm, hard cap 1 µm and at most one tenth of the shortest neighbouring edge or narrowest face. A sample or midpoint is not a proof.

Both candidates must reach the existing boolean, topology and export paths through the unchanged verification CLI. Incomplete requests, panics and missing certificates are not geometric refusals or successes. Uncertain decisions cannot be reported as proven.

## Fixed cases and baseline

Baseline is `main` at `d10c5c577d35cf3ae721275baee88a1017a634fa`, built before changing code. Preserve the executable SHA-256, doctor output, input hashes, complete responses and STEP files outside the repository. Freeze the ordered input list before candidate runs.

P0 coverage includes equal-diameter cylinders crossing at right angles, cylinders tangent along a line, a cylinder and sphere with equal radii, and poses along the cylinder seam. Exercise union, subtraction and intersection, including both operand orders and translated or rotated placements. Replay every equal-radius refusal and the cylinder seam requests from the existing hard-case survey, without publishing private inputs. Tangent non-manifold outcomes remain separate from numeric coverage failures. Replay the earlier-answer corpus and failure packs as regression evidence.

## Measurement and decision

Use `.claude/skills/verify-catallax/scripts/control-catallax.py` with `catallax-verify`. Do not change the CLI, its work allowance or its driver. Record fresh per-case elapsed milliseconds, refusal code and latency, valid topology, STEP presence, and certified deviation from the retained definition. An exact unchanged definition has zero geometric deviation; export deviation is reported separately. Sampled field agreement is not a global deviation certificate.

Compare manifests and strict STEP bytes for earlier proven answers. Check independent runs at one, two and four workers. Use pinned alternating matched pairs, one warmup and three measured pairs, for every changed outcome or timing claim. Keep machine contention with the results. The boolean click budget is 500 ms; recipe builds have 1,000 ms and complete comparisons 2,000 ms. An unavailable application loop is explicitly unmeasured.

Promote the candidate answering the most fixed cases inside the click budget with valid topology and no changed proven answer. A tie goes to implicit because it approximates nothing. Neither candidate qualifies without its certificates and regression checks. Verdict is **promote**, **park**, or **repeat once** only if the measurement was wrong. Promotion is a separate proper implementation on `main`.

## Results

**Verdict: park both.** Neither candidate established an inside-budget improvement with complete earlier-answer preservation. The implicit tie-breaker cannot waive those conditions. No repeat is requested; censored requests are recorded below, not mistaken for geometric refusals.

Both modes were built at `prototype/hard-intersections` head `6e4b7ab628cfcea71265ef075a78e27acdab7922`, using compile-time `CATALLAX_HARD_INTERSECTIONS=implicit` or `certified`. That branch is local-only, never merged and never proposed for merge. The verification CLI and driver are unchanged from the named baseline in both executables.

The implicit implementation evaluates plane/conic meetings directly in their exact frame instead of expanding a three-variable elimination polynomial. The certified implementation retains each regular closed chart beside chord tubes with a proven same-parameter error bound of at most 100 nm; it uses those tubes only to prove exclusion. Ambiguous tubes and singular branch charts use the existing exact solver. Neither implementation publishes approximated model geometry or establishes coverage for general singular meeting curves.

### P0 round

The frozen set contains 24 canonical cases across the four families above, all three booleans and both operand orders, plus 32 survey requests. Input SHA-256: `da26463e327095d94c695055db2bf1731c85a0077090d23b581adcf1eddd9530`. Runs used the default 100,000,000 work allowance, nice 10 and pinned cores on a shared Linux x86-64 host. External guards were 30 seconds per P0 request and 5 seconds per earlier-answer request. A guard expiration means incomplete, never refused.

| P0 metric | Baseline | Implicit | Certified |
| --- | ---: | ---: | ---: |
| Answered | 15 | 10 | 7 |
| Typed refusals | 3 | 2 | 0 |
| Incomplete at external guard | 38 | 44 | 49 |
| Answered latency median, ms | 969.493 | 857.719 | 11,126.962 |
| Answered latency maximum, ms | 24,533.600 | 25,896.379 | 24,960.966 |
| Typed refusal latency range, ms | 2,462.041 to 21,975.707 | 14,056.954 to 24,968.027 | Not returned |
| Valid exported STEP / returned STEP | 13 / 13 | 8 / 8 | 5 / 5 |
| Published model deviation, nm | 0 | 0 | 0 |
| Maximum reported STEP deviation, nm | 0.07484 | 0.01624 | 8.34e-9 |

Baseline refusals were `numbers.operation-budget` and `boolean.self-touching-shell`; implicit returned only the latter. Fewer returned refusals did not establish better coverage because incomplete requests increased. Empty answers correctly had no STEP. Every returned STEP passed OpenCASCADE validity and matched the manifest's face count. This checks exported topology, not unobserved native-body invariants. All returned models remained exact; the export bounds concern different completed subsets and are not an accuracy ranking.

Single-pass latency medians exclude incomplete requests and do not establish relative speed. The following pinned AB/BA pairs each had one warmup and three measured repetitions. One-minute host load ranged from 17.32 to 20.12 during these pairs.

| Matched P0 evidence | Implicit | Certified |
| --- | ---: | ---: |
| Fully paired cases with identical manifest and STEP | 11 | 7 |
| Median of per-case paired latency deltas, candidate minus baseline, ms | -90.046 | +13,558.748 |
| Candidate cases with paired median at most 500 ms | 0 | 0 |

The implicit mode completed one additional equal-radius cylinder/sphere subtraction during its single pass. Its separate pair attempt retained a warmup and two measured pairs with identical answers and STEP, but hit the 300-second external guard before the third pair completed. Both measured answers took over 39 seconds. It supplies no click-budget win or complete timing comparison.

### Earlier answers and limits

The earlier corpus and failure packs contributed 477 baseline attempts: 177 answered, 52 refused and 248 incomplete. All 229 baseline-complete cases were selected before candidate replay. Implicit preserved the complete manifest and strict STEP of 148 earlier answers, with 29 incomplete; certified preserved 70, with 107 incomplete. No completed earlier answer differed. This is not a complete no-regression proof. On 12 independent replay probes, each mode produced identical manifests, strict STEP and deterministic work at one, two and four workers.

The application loop was not exercised. The general and singular-branch cases still need a different construction or certificate before either idea can be promoted. Keep these implementations parked; any promoted replacement must first close the incomplete-answer and topology evidence gaps.

Executable SHA-256 values, with full inputs, commands, responses, pair samples and STEP artifacts retained outside the repository:

- Baseline: `ba5e92f873dedab5276bab82b1d0eb44db9f0cd36ac142711247e9d67e0cb136`.
- Implicit: `2a5e1a1b119d0285c1acdc54df9c84030715b39d1e66ea5d812df1008599aace`.
- Certified: `1550d72bfae8ed7e1430b48a38a08c5783fc1100ecec0588ee2ad000346eed32`.
