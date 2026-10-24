# Piece 1 contract tests

These tests were written from `docs/pieces/01-numbers-and-decisions.md`, the design and architecture documents and the number interface, without inspecting an implementation or reference kernel. Integer division follows floor semantics; case 7 and the interface-change log record the negative-divisor rule.

## Commands

CI executes all 35 blind contract tests, alongside the implementation's unit tests, the independent oracle checks and the follow-up regression tests. Release mode keeps the four-million-question orientation family practical:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release
```

`blind_test!` enables every case. The three comparison 01 exclusions are restored after correcting the statistics test and optimizing dyadic rational arithmetic. The full constants/identities test completed in 50.90 seconds and Chebyshev degrees 1 through 64 with per-root refinement in 138.76 seconds in focused release runs. The comparison report records the baseline observation limits and operation-level before/after timings. Run the entire piece 1 suite with:

```sh
cargo test -p catallax-numbers --tests --release --no-fail-fast -- --include-ignored --nocapture
```

The zero-proof follow-up adds exact arithmetic in a real quadratic field, zero identities by polynomial reduction in one algebraic generator, and product-formula separation bounds from structural Weil heights. Generator reduction never treats a nonzero remainder modulo a potentially reducible defining polynomial as evidence of nonzero value. Inverses require a constant polynomial GCD, so they are valid at every root. Unsupported expressions retain the annihilator fallback. Every proof cache is local to one attempt within a decision; the default policy and public signatures are unchanged.

The added regressions cover quadratic root-branch selection, positive square roots, equivalent radicands, independent quartic witnesses with tiny nonzero perturbations, and reducible defining polynomials with a divisor that vanishes only at another root. Existing repeated-decision and eight-thread checks still cover deterministic outcomes and work.

Case 39 also covers a four-radical expansion, a square root already in a quadratic field, and the negative sign of a three-radical sum. The four-radical zero is checked under another square root to exercise definedness proofs, and tiny perturbations must keep their signs. The imprint regression combines `sqrt(5/9)` with both roots of `648x^2 + 1116x + 145`. Further tower regressions cover reciprocal norms, positive square-root selection, fourth roots, and a tower too large for the exact attempt that must leave enough budget for refinement. The tower attempt has at most eight generators, 100,000 recursive steps, and one eighth of the remaining decision work capped at 1,000,000 work units; these are attempt limits, not changes to `BudgetPolicy`.

Constructible quartic leaves are also embedded in the tower. Removing the cubic term gives either a quadratic in the square or a cubic resolvent with a positive rational root, which splits the quartic over a real quadratic extension. Each candidate is checked against the defining polynomial and selected through its exact isolating interval. All four real roots of the irreducible quartic from edge-hunt case 0000163 are compared with independently constructed nested radicals, including tiny positive and negative perturbations. A second irreducible quartic has a nonzero linear term after depression and exercises the resolvent path. Unsupported quartics retain the existing exact fallbacks.

Working-precision enclosure tests include nested zero square roots and a positive `10^-100` radicand whose coarse ball overlaps zero. Bounded exact tower proofs preserve zero radicands without treating a small positive value as zero; unsupported overlapping radicands retain nonnegative-part interval square roots. A separate smoke probe records the precision loss before this fix and exact zero enclosures afterward, with unchanged decision statistics.

The working-precision suite also encloses `1/3`, `sqrt(2)` and the real root of `x^3 - 2` at 32, 64, 128 and 256 bits. Exact endpoint inequalities certify containment and widths must shrink as precision increases. Enclosures and expired requests preserve pre-existing decision statistics.

`root_isolation.rs` covers case 40's scaled large-leading-coefficient quartic and the exact probe-plane quartic from piece 5 case 23 at scale `10^-30`. Both must return four simple irrational roots under the default budget. Independent GMP Sturm certificates check the intervals. Further cases retain exact rational roots when denominator divisors cannot be enumerated, separate overlapping intervals across square-free factors, exercise the zero-precision bisection fallback, and refine roots at scales `10^-30`, `1` and `10^30`.

Independent directed-rounding and Sturm-oracle checks run in default CI. To run the directed-rounding check separately:

```sh
cargo test -p catallax-numbers --test enclosures oracle_directed_rounding_smoke
```

The benchmark retains its explicit opt-in:

```sh
CATALLAX_NUMBERS_RUN_BENCH=1 cargo bench -p catallax-numbers --bench speed
```

It prints decision counts, filter share, nanoseconds per decision by proof path, maximum precision and work, root timings through degree 100, three-sphere certification timings, pi timings at 1,000 and 10,000 decimal digits, and queries over 100,000 facts. There are no timing or filter-share thresholds. Refusals are reported, not hidden.

## Mathematical oracles

`rug` is a dev-dependency only. GMP rationals and integers compute exact reference arithmetic independently of `catallax-numbers`. MPFR directed rounding computes certified reference intervals for elementary functions. Oracle inputs come from fixed seeds and exact rationals, never from an implementation's answer. Tests compare exported exact numerators and denominators against the oracle.

The test-local Sturm sequence checks that an irrational root's defining polynomial is primitive and square-free, neither interval endpoint is a root, and exactly one root lies between them. Its arithmetic uses GMP rationals, not `UPoly`.

The elementary-function oracle uses greater precision than the requested result and requires the returned ball to contain the reference interval. The pi oracle uses 12,000 bits, exceeding 3,000 decimal places; e and ln 2 use 4,096 bits, exceeding 1,000 decimal places. These are independently computed certified constants rather than rounded floating-point literals. Identities with exact results also check zero, one, minus one and two directly.

The interval-sign regression checks addition, subtraction and multiplication for positive, negative and zero-crossing hulls and signed point intervals at 1, 3 and 128 bits. Midpoint-radius arithmetic, sign-stable products and direct dyadic rounding avoid redundant rational operations without changing rounding, including negative ties toward positive infinity. A separate before/after smoke probe compared 18,141 midpoint/radius/precision rows at 1, 2, 8, 64, 128 and 512 bits; every row was byte-identical. The full default-deadline answer for edge-hunt case 0000088 also retained byte-identical property enclosures from the earlier 300-second diagnostic run.

Dyadic sums/differences normalize the integer numerator before constructing the final denominator; dyadic products cancel powers of two before multiplication. The GMP regression checks canonical numerators and denominators at 1, 53, 2,048 and 10,000 bits, with negative denominators, zero, scaling and cancellation. Non-dyadic fractions retain the general GCD path.

Shared odd-denominator regressions check addition, subtraction, cancellation to zero and reciprocals against GMP at 53, 607 and 2,048 bits. Addition cancels only against the shared denominator GCD, and a reciprocal swaps already coprime fields without another GCD. Sign-stable ball products form their exact center and radius before rounding; a before/after arithmetic probe retained the same midpoint/radius hash across 40,000 products.

Equal non-dyadic denominators need only the numerator sum or difference and its final normalization, without denominator GCDs or cross-products. Zero and unit rational products reuse their normalized answers before cross-cancellation. The existing GMP field-law, shared-denominator, signed-cancellation and zero-domain tests cover these paths.

Real-expression construction evaluates operand-derived coarse balls only when no exact rational cache exists. The rational cache already takes precedence, so computing both balls discarded work. Expression dependencies, approximate intervals and definedness flags are unchanged. Existing shared-zero and undefined-outer-operation regressions cover those dependencies.

Ball construction rounds `n/d` with exact binary quotient arithmetic, then rounds the exact residual `|n - mid*d|/d` outward for its radius. Integers use the same binary centering path as arithmetic. Neither path needs to normalize temporary rational fractions. `ball_construction_preserves_signed_ties_and_exact_residuals` checks zero, positive and negative midpoint ties, precision-zero clamping and non-dyadic residuals; the existing GMP and MPFR families cover signs, scales and deep precision.

Dyadic alignment borrows the unshifted mantissa and shifts only the operand with the larger exponent. Equal exponents need no shifted copies. An unchanged rounded midpoint has exactly zero rounding error, without subtracting equal mantissas.

Ball products form the same exact interval before rounding. For sign-stable midpoint-radius inputs `(m,r)` and `(n,s)`, its center is `mn + sign(mn)rs` and radius is `|m|s + |n|r`. When only one interval crosses zero, the stable interval's endpoint farthest from zero scales the crossing interval's center and radius. Two zero-crossing intervals retain the four extremal products. These identities avoid expanding wide endpoints and then cancelling them again; midpoint rounding, outward radius rounding and work charges are unchanged.

`ball_alignment_handles_signed_exponent_gaps_and_exact_cancellation` compares signed, scaled addition, both subtraction orders and multiplication with GMP at 1, 53 and 1,024 bits. Exponents range from -257 through 257, including equal exponents, gaps in both directions and exact cancellation. The existing interval-sign family covers positive, negative, zero-crossing and singleton intervals.

Generated expression trees have depth at most six. Each carries an independent exact rational answer or proof of undefinedness. Square-root branches use `sqrt(v*v) = abs(v)`, or deliberately negative arguments. This supplies proofs of exact zero and undefinedness that a finite-precision oracle cannot supply. Separate irrational cases cover hidden identities, near-ties, algebraic roots and system-root coordinates.

All random families use the test-local SplitMix64 generator with a fixed seed per family. Numeric oracle answers do not depend on the clock. In-progress cancellation regressions run in debug and release builds without elapsed-time assertions. The resultant test uses its existing callback to set an already-expired deadline after the 32 Sylvester rows and the start of Bareiss elimination. Root isolation, pi and arctangent use a thread-local test clock. It stays just before the deadline until several checks in the targeted inner routine have succeeded, then advances to the deadline. Each call must return only `numbers.deadline`, never a partial answer.

The optional `test-deadline` feature exists only for process tests and is enabled by `catallax-cli`'s dev-dependency, not its normal dependency. With this feature, `CATALLAX_TEST_DEADLINE_CHECKS=N` selects the same clock and expires it on numeric work check N. The CLI regression uses check eight for all 18 sphere/cylinder boolean cases and verifies refusal-only manifests and removal of stale STEP output. Its 60-second process watchdog detects hangs, not deadline latency. Without the feature the environment variable is ignored, and all deadline reads use `Instant::now()`. Do not enable this feature in production builds.

## Coverage

| Contract cases | Tests |
| --- | --- |
| 1–7 | `exact.rs`: grammar and error precedence, 256-bit field laws, 10,000-digit integers, floor division with both divisor signs, binary ties/subnormals, decimal bounds |
| 8–13 | `enclosures.rs`: adjacent-double enclosures, outward arithmetic, overflow, underflow, zero and domain edges |
| 14–19 | `enclosures.rs`: directed high-precision constants and elementary functions, identities, precision changes, deterministic midpoint/radius |
| 20–31 | `decisions.rs` and `polynomials.rs`: complete 2D grid, translated near-planar 3D determinants, concave square-root signs, hidden zeros and algebraic equality |
| 32–36 | `decisions.rs`: tiny/random budgets, certified uncertainty bounds, depth-six expressions, all deadline-bearing entry points, repeated and eight-thread decisions |
| 37–39 | `decisions.rs`: hidden zero divisors, negative square roots, zero square roots and undefined-subexpression postorder |
| 40–46 | `polynomials.rs` and `root_isolation.rs`: exact rational roots, multiplicities, irrational isolation, scaled large-coefficient quartics, Wilkinson, Mignotte, Chebyshev and generated factored polynomials |
| 47–51 | `polynomials.rs`: one/no/multiple/singular roots, closed boundaries, both three-sphere intersections |
| 52–56 | `records.rs`: closure, symmetry, directed incidence, origin-insensitive deduplication, order and 100,000 facts |
| 57 | `records.rs`: mark order, worst bound by kind, exact/marked answers and outcome mapping |
| 58 | `decisions.rs` and `polynomials.rs`: every integer scale exponent from -30 through 30 |
| In-progress deadlines | `src/tests.rs`: root isolation, exact resultant elimination, high-precision pi and arctangent series; no partial proof or enclosure is returned |
| Shared expression definedness | `src/tests.rs`: a 64-level shared rational zero is proven, while a zero product retaining an undefined divisor still refuses |

The 2D family alone asks 65,536 sign questions at each of 61 scales, or 3,997,696 questions. The interval family adds at least 24,576 exact containment questions. The minimum of 30,000 generated questions therefore does not depend on counting helper assertions, repetitions, or benchmarks.

Additional polynomial tests cover canonical order, derivatives, substitution, resultants and square-free factors. Non-obvious sign and root expectations have derivations beside the cases. Tests do not prescribe a zero-proof strategy except where the contract explicitly requires exact arithmetic, and do not prescribe free-text refusal details.

Root-isolation regressions also cover rational roots `1/3` and `2/3` under 10,000 work units, rationally translated irreducible quartics, and the captured rotated cylinder-probe quartic. Rational polynomial signs try bounded exact evaluation before escalating interval precision. Rootless finite-field reductions exclude rational cubic roots of the quartic resolvent without isolating that auxiliary polynomial; primes dividing a coefficient denominator are skipped. The captured quartic previously exhausted 999,949,213 work units; the corrected run proves all four roots with 3,816,193. Independent GMP Sturm counts still verify every interval.

Low-degree isolation uses exact linear division and quadratic discriminants before generic square-free decomposition. A nonsquare discriminant excludes rational roots without denominator-bound reconstruction. Quadratics retain their earlier coarse brackets when bounded divisor enumeration was complete. Otherwise a stable quadratic formula avoids cancellation, derives the second root from the root product and certifies both branches with endpoint signs. Both endpoints must produce the same guard-precision midpoint/radius ball before the bracket is published. A short, work-bounded Descartes exclusion checks rootless quartics before a square-free GCD or auxiliary cubic; an inconclusive exclusion leaves the ordinary isolator responsible.

Guarded brackets do not replace an already tighter rational-reconstruction target bracket: its published enclosure remains unchanged. Strict signs can instead use the certified open root interval, including a zero endpoint. Positive and negative tiny cubic roots exercise that proof with no precision refinement and one work unit; a root at zero must still report coincidence.

The added regressions require scaled linear roots under 1,000 work units and large-denominator rational and irrational quadratics under one million. They check negative leading coefficients, distinct-root order, repeated-root multiplicity and nonreal roots. A product of positive quadratics must prove rootlessness under one million, while a tightly clustered real pair retains its quadratic factor and independent GMP certificates. Prime-only trial division keeps the original cutoff of 1,000 and the terminal prime certificate below 1,000,000; the boundary regression contrasts `999983` with `1009²`. An ill-scaled quadratic checks both positive roots without losing the tiny branch to cancellation. The CLI's public fourteen-step SolidDraft slot-cut fixture exercises the full subtraction, not only the standalone slot extrusion, under the unchanged 100,000,000 request allowance.

The translated clustered quadratic regression uses `(x-h)^2 - 2s^2`, with `h = 434139/7 + 1/(2^767-1)` and `s = 1e-20`. Its two roots are `h ± s sqrt(2)`. The original numeric baseline refused at 999,914,807 work units; the first exact-discriminant closed-form repair proved both roots with 900,624. The analogous cubic `(x-h)^3 - 2s^3` changed from refusal at 999,941,208 to a proven irrational root with 214,834,275. The rational cubic `(x-h)^3 - s^3` must still recover exactly `h+s`. GMP Sturm and normalized-power inequalities check the certificates without prescribing unrequested initial endpoint tightness.

Root-certificate reuse queries 257 distinct irreducible quadratics and revisits them in reverse order, exceeding the bounded memo's slot count. Refined signed-square inequalities and GMP certificates check each polynomial's own two roots across reuse and eviction. A separate warm-cache query preserves the original deterministic work charge and still refuses an expired deadline. Full-polynomial equality, not a hash alone, is required for reuse.

The certified-root enclosure regression isolates `sqrt(2)/(2^607-1)`, then requests width `1e-220` using a separate 10,000-work decider. A certified bracket already narrower than the working target needs no new endpoint-polynomial evaluation. The regression failed before that shortcut by retaining only a wider coarse enclosure; GMP signed-square inequalities check the narrow result independently. Calls that do need refinement retain their original bisection sequence and rounding.

Shared root snapshots memoize the original enclosure separately at each working precision. A cubic-root regression warms a certificate using an unrestricted working-precision call, then verifies that a limited-budget call still returns its original coarse ball and identical work statistics. GMP cubic inequalities check the warmed certificate; an expired-deadline request must still refuse. New intervals and changed defining polynomials start fresh memos.

The bisection-prefix regression varies precision history through 4,096 bits and work allowances from zero to one million. Cold and warmed cubic-root queries must return identical rounded enclosures, refusal codes and decision statistics, including when a cached prefix costs more than the remaining allowance. Independent GMP cubic inequalities check containment even when a zero-work query retains a wide coarse ball. A higher-precision request may reuse only an earlier prefix with the same sign-filter precision, never a later, tighter bracket.

Definedness checks reuse the immutable `fast_defined` construction proof instead of traversing shared expression graphs as trees. The regression proves a 64-level shared rational zero while retaining refusal for a zero product with an undefined divisor. Domain-error and in-progress root-isolation coverage remains.

`high_denominator_root_brackets_refine_with_a_small_work_budget` uses the strictly increasing polynomial `x^3 - 2` with valid endpoints `1 + 1/(2^4096-1)` and `2 - 1/(2^4096-1)`. It requires a `1e-30` enclosure under one million work units and checks its endpoints with independent GMP arithmetic. The unfiltered bisection path retained only a wider coarse ball. Budgeted interval signs avoid huge exact Horner products unless the sign remains unresolved.

The filtered-prefix cache regression warms the same odd-denominator cubic at 128, 512 and 768 bits, then compares cold and warm 768-bit policies under one million work units and the default work allowance. It checks full enclosure results and statistics across both the fixed-filter boundary and an underfunded cached prefix.

The content-normalization smoke probe uses a degree-16 monic polynomial with 65,536-bit coefficients. A unit leading coefficient proves the content is one independently of the implementation; all coefficients are retained. Existing `polynomials.rs` and unit tests check negative leading signs, rational content, common factors, zero polynomials and multiplicities.

`gcd_preserves_rational_factors_at_large_integer_scales` checks exact common factors when a finite-field image loses degree or cannot invert a rational denominator. A coprime pair with identical modular images must still have rational GCD one. These cases constrain the one-sided coprimality certificate without accepting a finite-field factor as a rational answer.

`singleton_zero_certificates_do_not_require_elimination` reproduces a discarded certificate using `0 * sqrt(2)`. One work unit must prove the known-defined zero, while an empty allowance still stops. A tiny nonzero product may be positive or unresolved, never coincident. A zero with an initially unresolved valid divisor must become provable after definedness is established; a hidden division by zero must still refuse despite its outer zero singleton.

Reciprocal balls use the exact center `m/(m^2-r^2)` and radius `r/(m^2-r^2)`, then the same nearest-midpoint and outward-radius rules, including midpoint rounding error. Integer quotient rounding avoids rational GCDs. The before/after smoke workload covers 4,806 signed, scaled and low-precision hull cases, including zero-crossing divisors and very small radii; all 40,000 results retain the same midpoint/radius fingerprint. The GMP/MPFR enclosure tests check containment independently.

Consuming temporary dyadic magnitudes leaves all arithmetic and rounding unchanged. The complete numbers suite covers signed products, reciprocal residuals and midpoint-rounding errors; a separate public jet smoke comparison checks 57,456 exact midpoint/radius/precision triples across orders 0, 1, 4, 8 and 20 at 8, 64, 128 and 256 bits.

The per-evaluation enclosure memo shares immutable balls internally and returns an owned root after releasing the memo. Existing clipping, undefined-node, nested-radical, root-cache budget-boundary and independent-oracle checks exercise local replacement without modifying another cached value. No shared reference escapes a numerical entry point.

### Pinned system coordinates

`system_certification_with_pinned_coordinate` certifies
`[x^2 + y^2 + z^2 - 1, x - y, z]` under the default budget, with `x,y` in
`[1/2, 1]`. It checks `z` in `[-1/2, 1/2]`, on either closed face, and in
the singleton `[0, 0]`. The returned coordinates must equal
`[sqrt(1/2), sqrt(1/2), 0]`, with independently checked enclosures.
`system_certification_with_two_pinned_coordinates` checks nonzero pinned
coordinates, including faces and a box with two zero-width axes. It also
checks exclusion and refusal for a box containing both roots.

Previously, intersecting with the Krawczyk image collapsed an exact linear
coordinate. Strict inclusion could never hold on that axis afterward.
Certification now requires closed containment of the image and the same
strict contraction bound as before. The Newton map is then a contraction
from the closed box into itself, so Banach's theorem proves existence and
uniqueness, including on a face or a zero-width axis. The full Jacobian
norm bound still proves nonsingularity. Exclusion, public signatures,
budgets and the piece 1 contract are unchanged.

## Verification limits

Comparison 01 observed 32 passing blind contract tests, two passing independent oracle checks and 16 passing implementation unit tests. One blind test failed an assertion and two did not complete within their observation windows. The follow-up corrects the test's undrained statistics call, adds deep-precision enclosure and GMP dyadic-normalization regressions, and restores both formerly unfinished contracts to default CI. All public signatures, proof methods and deterministic work charges are unchanged. A single local run does not prove cross-machine determinism.
