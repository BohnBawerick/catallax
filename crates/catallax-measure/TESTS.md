# Measurement blind tests

Source contract: `docs/pieces/06-measurement.md`. The suite was written blind using the documented interfaces and pieces 1, 2 and 3. Expected answers do not come from measurement, intersection or boolean implementations, reference kernels or an external tester.

## Commands

Ten complete blind bodies now run by default, including all six bodies repaired after comparison 06. The boolean cross-check alone remains ignored because its Steinmetz fixture refuses `boolean.imprint-ambiguous` before measurement. The 13 original private tests, two analytic-arc regressions and both independent fixture tests remain active. No proven-outcome assertion was weakened.

Run the active tests and independent derivations:

```sh
cargo test -p catallax-measure --release
python3 crates/catallax-measure/tests/derive.py
cargo bench -p catallax-measure --bench blind
```

To explicitly run every test, including the generated families and boolean cross-checks:

```sh
cargo test -p catallax-measure --release -- --include-ignored --test-threads=1
```

The inclusive command still exposes the out-of-scope piece 5 prerequisite refusal. The original comparison used a 1,800-second total observation and separate 600-second observations for bodies not reached by that run. Those historical observation limits are not test assertions or numerical refusals.

Never use debug mode for the full workspace suite. Piece 1 is expensive in debug mode.

Measurement shares each edge's trace across incident faces, and `mass_properties` shares those traces across volume, area, centroid and second moments. The property work budgets remain independent; the six second-moment components divide the inertia allowance equally. Cached Taylor enclosures are mathematical inputs, not cached measurement answers; standalone-versus-grouped consistency and deterministic enclosure checks still apply.

The inertia checks cover cylinder, sphere, cone, frustum and torus closed forms, exact rotated and translated blocks with signed products of inertia, and a disjoint boolean union against the sum of both blocks with parallel-axis shifts. Empty bodies have no tensor. The limited-precision regression checks the tensor's maximum entry width as the fourth property mark. Trigonometric antiderivatives now cover all monomials through degree six for trimmed torus second moments.

Compiled Taylor programs skip exact-zero convolution products, reuse each coefficient denominator's reciprocal across its jet, and clone inputs without first allocating zero jets. Symmetric products are computed once, but every addition/subtraction retains its original order so rounded endpoints stay unchanged. The polynomial-jet regression checks squares and positive square roots with zero coefficients between nonzero terms. Panel order, requested accuracy and budget defaults are unchanged.

Symmetric Taylor panels evaluate the center series only through order 18. Center orders 19 and 20 do not contribute: the order-19 moment vanishes, and order 20 comes from the whole-panel remainder, not the center expansion. The full order-20 remainder and its order-21 curve jets remain unchanged. Torus lifts still evaluate their full angle series. Per-operation work charges are unchanged; omitted computations no longer consume work.

`symmetric_panel_retains_the_signed_order_twenty_remainder` checks the analytic integrals on `[-1,1]` when the center coefficients stop at order 18 and the order-20 remainder is positive or negative. This catches accidentally dropping that remainder while trimming the center series.

## Coverage

| Contract | Test coverage |
| --- | --- |
| Cases 1-6 | Both blocks, cylinder, whole sphere, frustum, apex cone, loopless torus and Empty; volume, area and centroid; exact endpoints where required |
| Cases 7-13 | Hand-built through hole, reversed void shell, hemisphere, sphere cap, octant, both half cylinders and four-face Steinmetz solid |
| Case 14 | Every face has its own closed-form oracle; interval sum contains body area and overlaps the independently checked body enclosure |
| Case 15 | No invalid solids or planted faults; fixture topology and endpoint incidence are checked before measurement runs |
| Case 16 | Relative accuracies 10^-6, 10^-12 and 10^-18 on all specified curved cases; 10^-24 and 10^-30 on sphere, frustum and positive half cylinder; absolute floor 10^-40 |
| Case 17 | 64-bit/10,000-work and 64-bit/zero-work policies; exact rational results remain exact; other accepted answers must enclose the closed form; every missed width has its precise mark |
| Case 18 | Every case 1-13 variant at each decimal scale from 10^-30 to 10^30, plus 20 fixed-seed rational poses |
| Case 19 | Fresh-Decider repeats, fresh models on eight threads, identical mass-property and face endpoints; standalone volume/area agree with mass properties |
| Case 20 | Expired deadlines for all four entry points; either the ordinary answer with identical marks or exactly `numbers.deadline` |
| Case 21 | All nine rows of the boolean cross-check table; case 41 checks only its specified volume |
| Other rules | Signed and zero-containing intervals; default accuracy exercised on a sphere; body's marks stay on the body; thin-shell centroid separation; actual width and maximum centroid width in mark order |

The scale/placement/determinism family completes all 7,007 mathematical enclosures, excluding its additional standalone and face checks. There are fifteen solid variants with five quantities and one Empty variant with two. Each runs at 61 scales, 20 poses, two serial measurements and eight thread-local measurements: `77 * 91 = 7007`. The random-pose seed is `0x91a760c4f033e5b1`. Standalone volume now refines past the absolute floor until separated or budget-limited, so it agrees with the public volume returned beside a Solid centroid.

Only three bodies use the two tightest accuracy settings in correctness tests. There is no time-based assertion.

## Oracles and construction

`tests/common/oracle.rs` uses dev-only GMP rationals and MPFR directed rounding through `rug`. It never uses `Ball`, surface evaluation, traces or measurement to produce expected values. Rational results compare exactly. For transcendental results it first checks overlap, then refines the reference interval until it lies wholly inside the returned enclosure and is less than 1/1024 of its width. This handles negative half-cylinder centroids and subtractive hole volumes without decimal tolerances. An oracle ceiling of 32,768 bits reports an inability to certify containment rather than accepting overlap alone.

`tests/derive.py` independently checks the non-obvious formulas using Python `Fraction`: cross-section integrals for frusta, spherical caps and octants; the half-disc first moment; Steinmetz square sections; void and hole corrections; weighted cube centroids; the cylinder/sphere intersection; and the two-hole inclusion-exclusion volume. Alternating-series Machin bounds and integer square-root bounds also check pi and the frustum's irrational area. It imports no kernel or geometry library.

The hand-built fixtures use `add_*`, explicit rational vertices and recorded construction incidence. The void reverses both face senses and boundary travel. The upper and lower half cylinders keep the same positive-x seam. The octant has quarter-circle boundaries through a sphere pole. The Steinmetz ellipses are split at `(0, +/-1, 0)` with all four cylindrical patches oriented outward. Ring boundaries and loopless faces stay ring boundaries and loopless faces after scaling and placement.

Every fixture checks `Model::check`. Ungated checks additionally verify endpoint incidence with independent rational curve and surface equations and assert the intended topology counts at identity, tiny and large posed scales. They also cross-check the sum of the face formulas against the body's formula.

The rule-4 exactness requirement is applied only to the explicit planar bodies. Boolean planar results may contain meeting vertices, so their tests require containment and accuracy, not collapsed intervals. No test requires nonzero width or nesting between successive requested accuracies.

## Refusals and marks

Accepted budget-dependent upstream refusals use named codes, not prefix matching: `numbers.budget-exhausted`, `numbers.precision-limit`, `numbers.divisor-not-separated`, `numbers.sqrt-not-separated` and `numbers.root-not-certified`. Piece 4 can propagate a `NotExact` code as a refusal when it cannot locate a trace endpoint. Invalid-constructor and unsupported-geometry codes are not accepted for these valid supported fixtures.

Zero-work measurement refusals must be `measure.budget-exhausted` or a permitted upstream code. Mass-property budget checks, including the thin spherical shells, also allow exactly `measure.centroid-not-separated` when the budget cannot separate the public volume from zero, as rule 5 requires. No other refusal code was added. A stable analytic formula may separate the volume even under a small budget. Every returned solid centroid requires a strictly positive volume lower bound, including when volume is far below the requested absolute floor.

`measure.not-solid` and `measure.unsupported-surface` are unreachable through today's public body and surface variants. The suite does not fabricate future enum variants, invalid handles or malformed solids. These are not missing executable cases.

All returned marks must be `measure.accuracy-not-reached`, kind `Value`, one per missed property in volume/area/centroid order. The scalar bound is its actual width. A centroid mark uses the largest of all three coordinate widths, even when only one coordinate missed. Marks from the body's construction must not appear in measurement results or enlarge their enclosures.

## Analytic arc regressions

Polynomial boundary forms on rationally framed analytic arcs now use exact antiderivatives rather than Taylor panels. Tests cover quarter turns and complementary three-quarter turns with shifted circle reference directions. A sphere-meridian regression checks that a zero numerator cannot erase a gauge denominator whose pole lies on the closed arc. Unsupported radicals, variable divisors without a proof, and alternate sphere gauges retain certified Taylor quadrature.

Taylor products and quotient recurrences bound their loops by the last coefficient not proven to be the singleton zero. Skipped terms previously performed no charged arithmetic, so coefficient order, outward rounding and measurement work counts are unchanged. Coefficient checks use `Ball::is_exact_zero` without copying the midpoint and radius. Dense series retain the same recurrence.

## Benchmarks

Opt into timings:

```sh
CATALLAX_MEASURE_RUN_BENCH=1 cargo bench -p catallax-measure --bench blind
CATALLAX_MEASURE_RUN_BENCH=1 CATALLAX_MEASURE_BOOLEAN_BENCH=1 cargo bench -p catallax-measure --bench blind
```

The independent benchmark is named `blind` to preserve the creator's existing `speed` benchmark unchanged.

The benchmark requests each primitive and hand-built body at relative accuracies 10^-6, 10^-12, 10^-18, 10^-24 and 10^-30. The boolean-enabled command requests all nine cross-check rows, but stops if a fixture prerequisite fails. Output records wall-clock nanoseconds, all Decider statistics, reached widths, marks and refusal codes for volume, area, mass properties and each face. Centroid and inertia time are included in `mass_properties`; the reported widths include the six independent tensor entries. Body construction and oracle evaluation are outside measurement timers. The comparison report distinguishes completed cells, timed-out observations and shared boolean fixture failures.

Measurement's own Ball work and maximum precision are unavailable until an observation API exists; the reported counters are Decider-only.

The benchmark does not equate policy limits with work used, report fabricated zero counters, or gate speed. Without the environment opt-in it prints the availability notice and exits without running measurements.
