# Hard intersections

Protocol recorded before prototype code. Decision due at C2, 19 November 2026. Time box: one checkpoint. This is an experiment, not a change to the kernel's representation.

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

Not measured yet. The baseline revision and protocol above precede all prototype code. The result update will name both candidate revisions and executable hashes, give P0 time/refusal/topology/deviation counts against the baseline, identify unchanged-answer evidence and unavailable checks, and record the verdict. Prototype code stays on `prototype/hard-intersections`, is never merged and never opened as a pull request.
