# Piece 4 blind tests

The suite is derived from `docs/pieces/04-intersections.md`, the shared design and architecture contracts, and the finished number and shape interfaces. It does not use an intersection implementation or a reference kernel to compute expected answers.

## Commands

Run default checks:

```sh
cargo test --release -p catallax-intersect --all-targets
python3 crates/catallax-intersect/tests/verify_math.py
```

Every blind contract body is enabled. The complete hit scale/placement family passed a release run over all prescribed decimal scales and twenty fixed-seed rational poses before its ignore gate was removed. Its independent tilted-cylinder witnesses, exact residuals, contact labels and hit ordering remain unchanged. Default marks are still failures.

The two fixture-only tests exercise fixture construction and residual checks without calling any `catallax-intersect` function. Rational residuals are checked exactly. Radical residuals use certified enclosures as a smoke check, not as an equality oracle. The independent Python script checks non-obvious fixture arithmetic with `Fraction` and exact quadratic radicals. It also checks the two length constants numerically with separate smooth rational-parameter integrands and 50-digit Decimal arithmetic. Those numerical length checks are not claimed as certified quadrature.

Use release mode. The original inclusive run reached its 600-second observation cap; separate original bodies then had 180-second caps. Follow-up diagnostics retain every body and assertion, including a six-body sequential measurement on one logical CPU with 600-second per-body caps. A timeout is an incomplete observation, not a numerical refusal. The runtime target remains about ten minutes on one core; no test asserts an elapsed-time limit. The complete hit sweep is now verified under its original default policy.

## Coverage

| Contract cases | Files and checks |
| --- | --- |
| 1-27 | `surface_pairs.rs`, `common/cases.rs`, `common/oracle.rs`: every surface-pair configuration, both argument orders, exact analytic fields, irrational witnesses and forms, singular points and shared handles, contacts, point definitions, complete construction fact sets, store growth and repeated-call independence |
| 28-33 | `curve_surface.rs`, `common/curves.rs`: line, circle, ellipse and branch hits; empty and contained answers; multiplicity contacts; parameter-origin ordering; isolated or closed point definitions; own-surface fact decisions and refusal rollback |
| 34-39 | `curve_questions.rs`, `query_properties.rs`: membership and branch exclusion, curve equality with alternate witnesses and frames, line and cyclic order identities, analytic and one-sided singular tangents, strict arc interiors without adding the returned point to the store |
| 40-43 | `traces.rs`: whole circles and arcs, line segments, the quartic loop and both starts of the upper Viviani branch; endpoint and join enclosures, orientation, witness coverage, implicit polynomial residuals, Taylor identities and validated interval length sums |
| 44-45 | All 61 decimal scales from `10^-30` to `10^30`, 20 fixed-seed rational poses, both surface argument orders, and 30,000 additional rationally scaled and posed configuration questions. Curve-surface fixtures receive the scale and pose sweeps separately. |
| 46 | Eight-thread replay of surface intersections, curve-surface intersections and every curve query, including forms, tangents, arc points and traces; results, new definitions and decision statistics are compared. |
| 47 | Default, 64-bit/10,000-work, zero-work and expired policies; exact surface configuration with unchanged statistics; facts after expiry; allowed numeric refusals; bounded uncertainty and marks; no mutation on refusal |

The generated surface family cycles every table configuration. Rational scale factors vary dimensions and near-boundary distances; fresh rational poses vary axes and translations. Additional fixtures approach tangent and coincident boundaries from both sides by `10^-30`, include skew and oblique non-unit axes, and distinguish an irrational axis-length ratio from two irrational lengths with a rational ratio.

The surface oracle evaluates the published table algebra directly. It never asks `supports`, `intersect_surfaces`, `curve_form` or another piece 4 function for an expected answer. GMP integer square roots distinguish rational-square fields from intersection definitions. Swapped arguments are evaluated independently, including PP1 origins, CP3 order, CC6 roles and SS3's coordinate-dependent `perp`.

Returned meeting points are compared with `GeomStore::same_point` in a cloned store, so adding an expected point cannot affect the subject's construction facts. Geometric residual checks use bare positions where facts would otherwise make membership tautological. Tests do not pin proof kinds.

`NotExact` has no bound-kind field in the public interface. Query tests check the specified bound semantics numerically: distance to a dividing plane, or the magnitude of the ordering quantity. Operation marks do expose `BoundKind`; tests check their codes and kinds, merge-distance coverage, and the absence of a curve-incidence fact when branch membership was undecided. They accept either retained candidate after a merge, never a newly averaged point. Expired operations may give the exact answer or `numbers.deadline`, not a marked answer.

Trace tests do not prescribe a parametrisation or a piece count. They check orders 4 and 16 at 64 and 1,000 bits, including the polynomial identity `f(C(t)) = 0` coefficient by coefficient. Length bounds integrate speed enclosures obtained from `taylor(t, 1)`. Internal joins expose only balls, not exact point handles; their continuity is checked by overlapping enclosures. Endpoint orientation uses exact supplied endpoint handles and the arriving/leaving tangent appropriate to that end.

Trace constants retain exact `Real` definitions and cache enclosures on demand through `Decider::enclose_at`. Trace construction does not eagerly enclose constants to `max_precision_bits`. The trace regression checks that constructing a segment with an irrational endpoint leaves its constant enclosures empty, then evaluates Taylor jets at 64 and 256 bits and checks that constants use no more than the requested precision plus the evaluator's 64 guard bits. The numeric suite checks working-precision enclosure of known square roots at zero, divisor retries at higher precision, unchanged decision accounting and in-progress deadlines. High-order interval jets skip exact-zero arithmetic and reuse each denominator's reciprocal throughout a recurrence; the existing order-16 implicit identities and length enclosures exercise these changes.

Taylor recurrences also omit exact-zero coefficient tails and reuse symmetric products in squares and square roots without reordering rounded sums. A public-API before/after smoke probe covered 57,456 coefficient enclosures for circles and quartic cylinder/cylinder and cylinder/sphere branches, orders 0, 1, 4, 8 and 20, and precisions 8, 64, 128 and 256. Every midpoint, radius and precision fingerprint remained identical. The independent trace polynomial and length tests remain the correctness checks.

Trace clones share a bounded FIFO cache of Taylor enclosures keyed by the exact parameter ball, precision and order. Constant enclosures reuse the nearest cached working precision at least as high as the request, with outward rounding to the requested precision. A regression first refuses an unresolved reciprocal at low precision, separates it at higher precision, then checks that later lower-precision evaluations reuse that enclosure instead of refusing again. Cache warming can tighten later enclosures but never changes their validity. Both the circle Taylor regression and the constant-cache regression verify that an expired decider still refuses a cached request.

Surface-pair plans and rational Groebner bases have separate 32-entry thread-local FIFO caches. Keys contain complete ordered surface definitions or ordered polynomial coefficients, never store handles or addresses. Cached plans contain no store handles, and each intersection call still adds its own definitions and facts. Root isolation and every numeric predicate still run through the caller's decider. The `surface_definitions_and_argument_order_survive_repeated_queries` regression checks exact oriented sphere sections across fresh stores, changed centres and radii, swapped arguments, repeated calls and eviction.

Cache hits share immutable plans and basis slices instead of copying their vectors and arbitrary-precision coefficients. Constructing new store entities still copies the definitions those entities own.

## Resolved specification correction

Case 27's original description reversed the inner-tangency near misses. For cylinder radius 3, sphere radius 2 and centre `(d, 0, 0)`,

```text
x² + y² = 9
(x-d)² + y² + z² = 4
z² = -5-d²+2dx, with -3 <= x <= 3
max z² = 4-(3-d)²
```

At `d = 1 + eps`, this maximum is `4eps-eps² > 0`, giving CS6. At `d = 1 - eps`, it is `-4eps-eps² < 0`, giving CS4. The fixtures follow the mathematics and the CS4/CS6 table. The exact derivation is retained in the fixture comment and Python script.

## Non-gating benchmarks

Opt in to timing measurements:

```sh
CATALLAX_INTERSECT_BENCH=1 cargo bench -p catallax-intersect --bench blind_speed
```

The benchmark prints elapsed nanoseconds, decision counts, fact and filter shares, maximum precision and deterministic work. It measures surface pairs at identity and random poses, curve-surface pairs including fact-settled containment, curve queries, traces, and Taylor orders 4 and 16 at 64 and 1,000 bits. Near-boundary fixtures are labelled in the output. No timing result can fail a test.

The original benchmark currently panics when it tries to expose an analytic circle's mathematically rational parameter origin through `Real::as_rational()`. Comparison 04 retains all measurements before that panic and measures the four remaining query groups in scratch with their exact fixture origins. The report distinguishes that benchmark assumption from an intersection failure; the production benchmark is unchanged.

Without `CATALLAX_INTERSECT_BENCH=1`, the benchmark prints its activation command and exits without exercising intersection operations. This target also remains inert during `cargo test --all-targets`.
