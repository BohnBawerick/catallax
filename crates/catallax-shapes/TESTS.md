# Piece 2 blind contract tests

Written from `docs/pieces/02-shapes.md`, the design and architecture contracts, the finished `catallax-numbers` API and tests, and the interface-only `catallax-shapes` crate. No implementation branch, other worktree, or reference kernel was consulted. Public signatures and implementation bodies are unchanged.

## Run the tests

All 34 blind contract bodies are enabled. Two filter-dependent assertions use
the smallest sufficient budget, `max_work = 1`; assertions expecting
`NotExact` retain `max_work = 0`. The historical corrected run passed
those bodies, 13 implementation unit tests and the independent oracle.

Compile every test and benchmark, then run the enabled tests and independent
fixture-oracle check:

```sh
cargo check -p catallax-shapes --all-targets
cargo test -p catallax-shapes --lib --tests --release
```

During the implementation/comparison stage, run **all** contract tests, including the ignored tests, on one test thread:

```sh
cargo test -p catallax-shapes --tests --release --no-fail-fast -- --include-ignored --test-threads=1 --nocapture
```

The determinism test creates its own eight threads. There is no elapsed-time assertion. The corrected complete release run took 86.93 seconds across its test binaries, excluding compilation. This is one shared-host observation, not a speed threshold.

Check formatting and lints:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
```

## Oracles and generated questions

GMP rationals through the test-only `rug` dependency evaluate the implicit formulas and quaternion matrices independently. Fixed-seed SplitMix64 supplies inputs. Hand-answered cases compare exact polynomial coefficients, rational values, signs, or enclosure containment. Non-obvious answers have derivations beside the cases.

The main generated family asks exactly 32,940 geometric questions:

```text
20 rational poses × 61 scale exponents × (5 surface kinds × 3 offsets + 12 curve cases)
```

Every exponent from -30 through 30 occurs. The three offsets are on, inward and outward, with exact displacement factor `10^-30`. Surface points come from rational circle parametrisations. The expected sign comes from the table formula evaluated in GMP rationals, not `Surface::implicit` or any decision. Rigidly moved point coordinates are calculated by test-local rational matrix arithmetic.

Additional families check the named cases at every scale and at 20 poses, 64-bit rational quaternion identities, implicit invariance, parameter round trips, and repeated/eight-thread decisions. The larger sweeps use rational positions; the dedicated tests cover square roots and certified system roots without multiplying expensive hidden-zero questions across the full sweep.

The parameter liveness family cancels both circle and quaternion denominators when choosing its scale. Its posed points have integer coordinates, so the input balls really have radius zero. All chosen angles are in the first quadrant and separated from seams, axes and poles by more than the contract's `L / 100`. This avoids demanding liveness for a nonzero-width enclosure of a non-dyadic rational.

`oracle.rs` runs without an implementation. It checks the independent fixture arithmetic against the contract's quaternion matrix, Hamilton product, torus values, rational-circle identities and the strict rational bracket used by the uncertainty fixtures.

## Coverage

| Contract cases | Test files and checks |
| --- | --- |
| 1–10 | `poses.rs`: given matrices, quaternion scale/sign/hash equality, Hamilton products, orthogonality, determinant, inverse, composition order and point/vector translation |
| 11 | `analytic.rs`: every constructor refusal, direction/size/perpendicularity precedence, positive near-degenerate radii and axes, exact non-perpendicularity |
| 12–15 | `analytic.rs`: expanded coefficients without rescaling, posed implicit values, known moved cylinder, gradients, singular apex and every transformed field |
| 16–23 | `decisions.rs`: exact carriers, opposite cone nappe, near points and large translations, closed square roots, both three-sphere roots, tiny budgets, deadlines and undefined-coordinate refusal |
| 24–29 | `analytic.rs`, `decisions.rs`: equation order, analytic incidence, intersection carriers rather than branches, equal points represented differently, exact parameter origins and distance bounds |
| 30–36 | `decisions.rs`: equal surfaces written differently, near-unequal defining data, different kinds, normal signs and zero gradients, deadline/statistics exceptions |
| 37–44 | `parameters.rs`, `generated.rs`: values and derivatives, non-unit/reversed frames, advisory domains, output precision, interval inputs, seams and poles, off-surface inverse formulas, liveness and round trips |
| 45–49 | `store.rs`: handle order, no deduplication, typed key ordering, construction-fact order/origin, transitive fact paths, cached coordinates, exact positions and ordered enclosure calls |
| 50–52 | `generated.rs`: scale and placement families, all six decisions repeated, eight-thread equality of outcomes and statistics |
| Decision fallback rules | `uncertainty.rs`, `decisions.rs`: zero-work unresolved expressions, distance bounds for every surface and analytic curve kind, intersection-carrier bounds, normal-dot value bounds, nonzero/refusal short-circuit, unresolved code/detail, no fact mutation |

Tests prescribe a zero-proof kind only where the contract does: construction facts and exact-rational `same_surface`. A valid distance bound need not equal one chosen conversion formula. Undefined angular parameters may return `None`; seam tests constrain only defined values when `Some` is returned. No test compares approximate coordinates with an epsilon.

Point equality schedules cheap nonzero rational values and cached conservative interval signs first. When these do not separate, it requests metered enclosures at the existing interval widths before exact-zero comparisons. The chosen coordinate is still certified by the decider; an interval containing zero proves nothing. Existing late-coordinate separation, near-coincidence, zero-work distance-bound, undefined-coordinate and deadline cases remain. Tests do not prescribe the number of coordinate queries before a separation witness.

## Speed measurements

The blind benchmark is the `blind_speed` target. The creator's separate `speed` target is preserved unchanged. `blind_speed` exits before calling any implementation unless explicitly enabled:

```sh
CATALLAX_SHAPES_RUN_BENCH=1 cargo bench -p catallax-shapes --bench blind_speed
```

It prints nanoseconds per call, sample counts, outcome counts, filter/fact shares, maximum precision and work. Workloads cover:

- `side` by surface kind in general position, near points at `10^-30`, square-root positions and certified roots;
- fact-only and mixed decisions;
- `point_on_curve`, `same_point` and `normals_agree`, including intersection curves;
- `evaluate`, `derivatives` and `parameters_of` at 64 and 1,000 bits;
- `implicit`, transforms, quaternion construction, and point enclosures at widths `2^-64` and `2^-1000`.

Refusals and absent parameter results are printed. No timing or proof-path-share threshold gates a test.

## Verification boundary

The initial comparison's failed assertion was a test-budget error.
After correcting the two filter-dependent budgets, the body was enabled.
These are historical observations, not new publication-run results.

Case 28 passed with the merged piece 1 `certify_root` fix for a box where one
coordinate of the root is pinned exactly. Its specified box is unchanged.
