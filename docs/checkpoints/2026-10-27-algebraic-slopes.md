# Algebraic slopes diagnosis, 27 October 2026

This note diagnoses where 82° and 100° countersinks and drafted parts spend work. It measures only; it does not change any crate. All Catallax numbers below come from `main` at `e3a7f71ce05c742690a496353fc58cf1c89ec0e1`, built with `cargo build --release`. The `catallax-verify` binary had SHA-256 `de2375a67bb362e08d8802a4828185f762dee12dbad1400c98df24005fefada3`. The SolidDraft numbers come from SolidDraft `676b8e8357a09eeca936f4928ab5805fb21611ad`, which pins Catallax `a5c42583dd684e54c473bc377267e5e190c37bb6`. Work counts are deterministic for a build and input. Wall times were measured on a shared 4-core cloud container and are diagnostic only.

## Summary

| Finding | Status on `main` |
| --- | --- |
| The recorded topology query of about 4e8 on a drafted part | Not reproduced. The draft scenario's largest click-path call is 9.3e6. The scenario fits the default 1e8 allowance. |
| The recorded "more than 1e8 in the hole operation" for 82° and 100° | Not reproduced on a plate. The 82° hole costs 1.8e7 and the 100° hole 3.1e6. |
| The requested volume measurement of an 82° or 100° countersink | Over the allowance. It costs 1.1e8 and about 28 s. Cylinder ring integrals are skipped because of a syntactic rational gate. |
| A second 82° or 100° countersink in the P5 housing | Refused with `numbers.budget-exhausted` after about 1e9 work. Generic Gröbner coefficient zero tests cause this. |
| Strict STEP publication of a 3° or 5° drafted block in the verification exchange | Refused with `numbers.operation-budget` at 2e9. Decimal coordinate certification causes this. |

## Inputs

All inputs are public and synthetic.

- **Countersink plate.** SolidDraft's `ci/scenarios/countersink.json` builds a 60×40×5 plate with two through holes of diameter 4 at (15,20) and (45,20). Each hole has a countersink of diameter 8. The variants change only the countersink angle to 82° or 100°, at build and at both edits. The Catallax recipe that SolidDraft exported for the build step is the `catallax-recipe/4` question in that step's command bundle. Its hole step is `{"op":"hole","target":"op_3","depth":"through","diameter":"4","tip":"flat","head":{"countersink":{"angle":"82","diameter":"8"}}, ...}`.
- **Drafted block.** SolidDraft's `ci/scenarios/draft.json` builds a 40×30×20 block, drafts two sides by 3° about the top datum, adds a boss of diameter 10, edits the draft to 5° and undoes.
- **P5 housing.** This is `parts/holes/P5-countersinks.json` with the countersink angle set to 82° or 100°. Further variants keep one hole, two adjacent holes, two opposite holes, or two holes in separate hole steps.
- **P1 drafted.** These are `parts/P1-drafted-bracket.json` and, as a control, `parts/P1-bracket.json`.

The verification exchange was driven through `catallax-verify stream` with `source: "recipe"`. The SolidDraft click path was driven through `soliddraft mcp --max-work 100000000` with the scenario's own tool calls, `scenario_measure` and `scenario_evidence`. This container has no graphics adapter, so `scenario_prepare` and the viewport were skipped and no viewport timing is claimed. The measurement split by property used an untracked example that calls the public `catallax_measure` and `catallax_mesh` functions under `with_operation_budget`. That example was not committed.

## SolidDraft click path at the default allowance

Work per kernel call, at an allowance of 1e8 per call:

| Scenario step | Largest calls | Requested measurement |
| --- | --- | --- |
| Draft 3° build | draft 5.2e6, topology 1.3e6, preview 4.7e5 | mass properties 6.0e6 |
| Boss on the drafted block | boolean 9.3e6, topology 1.3e6 | mass properties 6.0e6 |
| Draft edit to 5° | boolean 9.0e6, draft 2.8e6, topology 1.0e6 | mass properties 1.8e6 |
| Countersink 90° build and edits | hole 1.6e6, preview 2.3e5 | completes, under 130 ms |
| Countersink 82° build and edits | hole 2.1e7 to 2.2e7, topology 2.2e6, preview 6.9e5 | **refused after 1e8 and about 28 s** |
| Countersink 100° build and edits | hole 4.6e6 to 4.7e6, topology 8.8e5 | **refused after 1e8 and about 27 s** |

The `soliddraft export` of the drafted block to STEP and STL took 0.7 s and 1.1 s, and `soliddraft draw` took 0.7 s. Neither CLI command reports work.

## Countersink measurement: work by property

These counts are for the 82° and 100° plates on `main`, with an allowance of 1e9 and one measurement per row:

| Property | 90° | 82° | 100° |
| --- | ---: | ---: | ---: |
| Volume | 332 | 2,690,864 | 1,353,878 |
| Volume, area and centroid | 1,156 | 4,713,648 | 3,397,166 |
| Mass properties with inertia | 5,073 | 109,855,772 (28.7 s) | 109,469,414 (25.9 s) |
| Topology face areas (SolidDraft `topology`) | 337 | 2,221,487 | 877,205 |

The face areas split as follows. At 82°, each cone face costs 2.8e3 to 6.4e4, while the cylinder faces cost 1.7e6 and 4.2e5. The planar faces cost under 60.

### Decision class

`catallax-measure/src/ring.rs` integrates a whole circle on a plane, cylinder or cone in closed form. Since #157, a cone admits its ring only after a decided coaxiality proof. The cylinder branches of `compile_admitted`, `second_moment` and `face_centroid` still admit a ring only if `(center - origin) × axis` has zero squared norm and the circle and cylinder have equal squared radii, both tested with `as_rational()?.is_zero()`.

The hole builder places the ring where the countersink cone meets the bore at the depth `(R - r)·tan(θ/2)`. The slope is `sin θ / (1 + cos θ)`, the cosine is a root of a degree-24 integer polynomial, and the sine is `sqrt(1 - c²)`. The ring's center coordinates are therefore expression trees such as `15 + 0·(...)`, and `as_rational()` returns `None` even though the offset is exactly zero. The gate then rejects the ring, and the cylinder face falls back to traced quadrature. Stack samples of the inertia run show this: all 25 samples are in `engine::integrate_face`, then `integrate::Program::panel`, then Ball arithmetic. At 90° every coordinate is rational and the closed form is used.

The decider proves the same offset zero (`ProvenCoincident`, algebraic proof) for 62k work at 82° and 3.1k work at 100°.

### Proposed exact reduction

The proposal is to admit cylinder rings with a decided zero test, as #157 did for cones. Where `as_rational()` gives `None`, the test calls `Decider::sign` on the squared offset and on the radius difference. It admits the ring only on `ProvenCoincident`, rejects it on `Proven`, and propagates `Refused`. A `NotExact` outcome keeps today's quadrature fallback. The integral is the same closed form, so the answer does not change. Only a proof that is already decidable is no longer discarded.

A local experiment applied this change to every remaining `as_rational()?.is_zero()` gate in `ring.rs`. All of those gates are cylinder gates on `main`. The experiment was then reverted. It measured:

| Property | 82° on `main` | 82° with decided gate | 100° on `main` | 100° with decided gate |
| --- | ---: | ---: | ---: | ---: |
| Volume | 2,690,864 | 77,975 | 1,353,878 | 19,309 |
| Volume, area and centroid | 4,713,648 (742 ms) | 108,141 (5 ms) | 3,397,166 (780 ms) | 49,985 (5 ms) |
| Mass properties with inertia | 109,855,772 (28.7 s) | 718,087 (89 ms) | 109,469,414 (25.9 s) | 330,515 (95 ms) |
| Topology face areas | 2,221,487 | 72,598 | 877,205 | 13,842 |

The 90° counts did not change. For 82° and 100°, every enclosure was compared: the volume, the area, three centroid coordinates and nine inertia entries. All 28 enclosures with the decided gate lie inside the enclosures on `main`, and they are four to six orders of magnitude narrower. The midpoints agree to every printed digit. For example, the 82° volume is 11797.237863199 mm³ and the 100° volume is 11818.099296802 mm³.

This reduction is in `catallax-measure`, outside the `catallax-numbers` and `catallax-boolean` cores. It needs no new public type. Its acceptance is the SolidDraft countersink scenario at 82° and 100° with the default allowance and its volume checks. Because enclosure widths shrink for these bodies, measurement-manifest preservation must be checked for every corpus row with an irrational cylinder ring. Rows whose rings are rational keep the syntactic fast path and are unchanged.

## P5 housing: a second countersink refuses

These are verification-exchange recipe requests with an allowance of 2e9, on `main`:

| Variant | Hole operation work | Outcome |
| --- | ---: | --- |
| One hole at 90° / 82° / 100° | 21,299 / 8,874,061 / 709,629 | answered |
| Four holes at 90° | 520,719 | answered |
| Four holes at 82° | 1,149,913,820 (76.6 s) | `numbers.budget-exhausted`, 2048 bits, 1 not exact |
| Four holes at 100° | 1,029,029,548 (0.5 s) | `numbers.budget-exhausted`, 16384 bits, 1 not exact |
| Two holes in separate steps at 82° | first 8,874,061; second 1,119,067,312 | second refused as above |
| Two holes in separate steps at 100° | first 709,629; second 1,008,205,373 | second refused as above |

The refusal comes from the decider's per-decision cap, not the request allowance; `exhausted` is false. All 15 stack samples of the 82° run follow the same path: `boolean::graph::intersect`, `intersect::hits::intersect`, `intersect::solve::intersections`, `algebraic_solve::intersections`, `groebner` and `reduce`, then `Polynomial::add`, `Decider::sign`, `separation_bound` and resultant-built annihilators. `solve::intersections` takes its rational path only when every equation is rational. The second hole's circles against the first hole's cone, both with coefficients in the tower generated by the angle root and its square root, go to the generic Gröbner basis. That basis tests every cancelled coefficient for zero with a separation bound over a deep expression tree. The faces are 24 mm or 34 mm apart.

There are two candidate exact reductions; both are in `catallax-intersect` and go through L01. The first rules out pairs whose curve range and trimmed face range are disjoint, using the existing rational box enclosures before any solve. The second keeps the coefficients in a single algebraic extension with a known minimal polynomial, so that a zero test is exact reduction rather than a separation bound. The first changes no answer, because it removes only pairs that have no admissible hit. It was not measured here.

## Drafted parts

The drafted block costs 4.0e6 for the draft step at 3° and 1.8e6 at 5°. At 45°, where the tangent is rational, the whole request costs 1.0e4. In the verification exchange, strict STEP publication of the 3° and 5° drafted block refuses with `numbers.operation-budget` at an allowance of 2e9. This happens with or without the boss. All samples follow `step::export::Exporter::curve_geometry`, then `number::point`, `number::enclosed`, `Decider::compare` and `Decider::sign`, ending in a separation bound. Here `enclosed` proves that an irrational coordinate equals the decimal written for it. SolidDraft's replay of its own exported draft question refuses the same way at the scenario's 1e9 allowance. The SolidDraft `export` command takes a bounded-deviation path and writes STEP in under a second. This path is in `catallax-step`, owned by L07, and this note proposes no change to it.

`parts/P1-drafted-bracket.json` costs 197,853,857 in total, and its undrafted control `parts/P1-bracket.json` costs 195,790,933. The difference is the draft, about 2e6. The rest is the slot Booleans: 3e6 to 1.6e7 each, and 1.1e8 for the final `slot-middle` subtraction. This cost does not come from irrational slopes.

## What remains

- The cylinder ring gate reduction: a `FULL SIX` slice in `catallax-measure` with the acceptance above.
- The SolidDraft draft scenario can drop its `max_work` of 1e9. That is a SolidDraft scenario change, and its owner should make it after confirming the full runner with a graphics adapter. The `S/docs/needs-catallax.md` lines on the 4e8 topology query and the 1e8 hole cost are stale for the pinned build. The measurement cost of 82° and 100° countersinks, recorded above, is the real gap.
- The P5 multi-countersink refusal and the strict STEP decimal certification of drafted coordinates belong to L01 and L07 respectively.
- These timings come from a shared container, so no speed claim follows. Matched pairs on a quiet pinned core must qualify any reduction.
