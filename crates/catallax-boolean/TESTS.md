# Piece 5 blind tests

The suite was written from `docs/pieces/05-boolean-engine.md`, the finished pieces 1 to 3, and the public piece 4 and 5 skeletons before either implementation was available. No boolean or intersection implementation supplied an expected answer. Post-comparison repairs preserve the public interfaces and mathematical result assertions.

## Run

Run the enabled real-build cases and compile the quarantined cases:

```sh
cargo test -p catallax-boolean --release
```

Ignore reasons distinguish unresolved implementation failures from opt-in long-run qualification. Mathematical expected answers are retained. The historical real-build inventory below records the original failures and observation timeouts.

Run every test, including ignored cases, the complete placement qualification, both refusal cases at nodes, and the eight-thread determinism checks. This audit also runs the deferred placement family and is not expected to pass until its implementation refusals are settled. `--no-fail-fast` ensures a failure in one test target does not skip the remaining targets:

```sh
cargo test -p catallax-boolean --release --no-fail-fast -- --include-ignored --test-threads=1
```

Run the deferred default placement family explicitly:

```sh
cargo test -p catallax-boolean --release --test properties case_57_default_scale_and_placement_variants -- --ignored --exact --nocapture
```

Case 57 prescribes three variants per operation: identity placement at exponents -30 and 30, then one fixed-seed exponent in -29 through 29 with pose `i mod 20`. The family contains 456 booleans. The 2026-10-20 shipment decision defers the entire default family to the next fix batch after three observed implementation refusals. All variants, counts, named points and chain checks remain in the ignored body; no refusal is accepted as a mathematical answer. The temporary per-cell exclusions and focused quarantine tests were removed.

The observed cells are `29-bite` Unite at scale `10^-30` and identity pose, and `23` Unite at scale `10^-20` under its seeded pose, both refusing at Imprint with `numbers.divisor-not-separated`; `23` Intersect at scale `1/10` under its seeded pose refuses at Sew with `numbers.budget-exhausted` (`root isolation exhausted its budget`). Their evidence files are `case29-scaled-divisor-after-quartic.json`, `case23-seeded-imprint-divisor-refusal.json` and `case23-seeded-intersection-sew-budget-refusal.json`. The latter run took 131.14 seconds including two preceding identity-placement cells, not as an isolated-cell timing or a full-default pass. Further placement investigation and divisor capture belong to the next fix batch.

A chain keeps the third variant of its first operation while `i` advances by every step. The ignored qualification test executes every exponent at identity and every one of the twenty poses at exponents -30, 0 and 30, for 121 variants per operation and 18,392 booleans, and prints its wall time:

```sh
cargo test -p catallax-boolean --release --test properties case_57_long_run_scale_and_placement_qualification -- --ignored --exact --nocapture
```

The generated material-identity test separately requires at least 30,000 checked containment questions against a few hundred boolean results. A counter increments only after an oracle-answered query is checked. Primitive boundaries are excluded from that generated result family; the named boundary queries run separately.

No test body asserts an elapsed-time limit. The historical real-build review ran each test in its own process with a 600-second observation cap. The cap bounded that review, not mathematical correctness. The full Cargo command above is uncapped. To replay capped execution, compile with `cargo test -p catallax-boolean --release --no-run --message-format=json`, enumerate test executables with `--list --format terse`, then run each name with `timeout 600s <executable> --exact <name> --include-ignored --test-threads=1 --nocapture`.

The periodic-chart regression builds a closed polygon from independently lifted segment boxes at 8, 64 and 128 bits. It checks winding after chart alignment. The propagated shift is the exact integer period, not the difference between rounded midpoints; the latter can make adjacent endpoint enclosures disjoint at every precision. Geometric direction searches have 64 attempts and inward probes have 256 dyadic depths, independently of `max_precision_bits`. Exhaustion retains the existing refusal paths.

## Coverage

| Description cases or rule | Tests | Assertions |
| --- | --- | --- |
| 1 to 9 | `contract.rs`, `common/cases.rs` | General, full and partial face contact, edge and corner contact, through/coaxial blocks, equal size, identical operands, both `10^-30` near-miss sides |
| 10 to 16 | `contract.rs` | Through and blind holes, external tangent line and both near misses, inscribed cylinder, seam groove, enclosed cylindrical void |
| 17 to 22 | `contract.rs` | Corner and face sphere cuts, outside tangency and both near misses, inside tangency, inscribed ball, seam and poles |
| 23 to 29 | `contract.rs` | Unequal crossing cylinders, Steinmetz ellipses, parallel overlap/tangency, both coaxial configurations, node contact, external point contact and tiny bite |
| 30 to 35 | `contract.rs` | Rational and irrational circle-plane heights, equatorial contact, Viviani node, one-branch bite, point contact, seam |
| 36 to 40 | `contract.rs` | Crossing, external tangency and both near misses, internal tangency, concentric and identical spheres |
| 41 to 44 | `properties.rs`, `common/chains.rs` | Every intermediate result's counts, topology, source surfaces, orientation, marks, reports, immutable predecessors, exact membership probes; two holes, four-block frame, two lens cuts, five unions |
| 45 | `contract.rs` | Both operand orders for empty identities, distinct empty handles, result kind, counts, marks, unchanged operands; either permitted report form |
| 46 | `contract.rs` | Meeting cone refusal `intersect.unsupported-pair` at `intersect`, both operand orders and all operations; a disjoint cone may refuse or return the mathematically correct result |
| 47 | `contract.rs` | Connected line-contact material refuses with `boolean.self-touching-shell` at `sew` |
| 48 to 50 | `stages.rs` | Two circle rings with exact centres/radii and owners; one isolated rational point; four square section edges with exact endpoints and owner lists in body order |
| 51 | `stages.rs` | Every pair variant's imprint pieces remain on an operand surface; exact face witnesses, independent containment and outward-normal classes, unrelated face returns `None`, common edges occur on both sides |
| 52 | `stages.rs` | Exact keep table, order and reversal senses for all three operations; direct sewing is checked against the same counts and material oracles |
| 53 | `stages.rs`, `common/mod.rs`, pair and chain runs | Fresh and seeded callers; five records in order; every additive counter preserved as `after = before + sum(records)`; nondecreasing running `max_bits` from the caller's prior peak to its final peak; mark counts and a decision-free keep stage |
| 54 | `properties.rs` | Primitive interior, exterior, face, edge, corner and rational rim queries, including the sphere's `10^-30` perturbations |
| 55 and 56 | `properties.rs` | Fixed-seed rational points, near-boundary probes, random primitive poses, at least 30,000 result queries, union/intersection commutativity and exact set identities |
| 57 | `properties.rs` | Default extreme-scale and seeded placement variants; opt-in 121-variant qualification; unchanged counts and transformed named points, including every chain intermediate |
| 58 | `properties.rs` | Two fresh runs and eight fresh threaded runs, equal counts, marks, statistics, handles, public topology and geometry records |
| 59 | `properties.rs`, `provenance.rs` | Prescribed small budget and expired deadlines across pair cases and chains; correct exact answers, sound marked answers, stable refusal codes and stages, immutable inputs |
| NotExact containment | `properties.rs` | Every observed NotExact retains `numbers.budget-exhausted` and has a nonnegative, mathematically valid Distance bound to a finite primitive face; all Proven outcomes remain correct |
| Marks and facts | `contract.rs`, `provenance.rs` | Stage marks in execution order, then first and second operand logs without deduplication; body and answer logs agree; required incidences have documented construction/proof provenance |
| Nonmaximal operands | `maximal.rs` | Equator-split sphere, retraced slit, one-vertex cylinder rims, subdivided block edge; removal of redundant faces, edges and vertices without changing operands |
| No rational point | `contract.rs` | Hand-built cylinder with no rational point refuses with `boolean.no-rational-point` |

Per the `boolean-stage-statistics` ruling, reports use `dec.stats()` snapshots without resetting the caller. All fields except `max_bits` are stage-local counter increases. `max_bits` is the running maximum at the stage's end, so keep has zero counter increases but retains the preceding stage's precision maximum. Case 53 runs every listed pair operation with a fresh caller and with prior exact-number decisions, including a narrow enclosure of `sqrt(2)`, to catch lost totals and incorrect stage-local peaks.

All listed pair results check lumps, shells, faces, edges, ring edges, vertices, surface-kind counts and curve-form counts. Curve kinds come from `catallax_intersect::curve_form`, not `EntityCounts::edges_by_kind`. This distinction is exercised by the irrational lines in case 13 and circle carriers in case 30.

`topology_rules` additionally checks that different lumps share no face, edge or vertex, that edges have two opposite uses in one shell, that no slit or same-surface face division remains, and that degree-two vertices between equivalent curves are removed. A singular vertex may occur twice within one face loop, but two loops of a face may not share it. Named node and contact counts check the required duplicated or absent vertices and ring edges.

`boolean.not-solid` has no constructible operand in the current `BodyKind`, exactly as case 46 states. Ambiguity refusals are checked when observed under restricted budgets; the suite does not require a particular budget to force a particular creator-chosen decision path. In particular, an exact rational method remains free to succeed at a small budget.

## Independent mathematics

`common/oracle.rs` evaluates box inequalities, finite-cylinder radial and axial inequalities, and sphere quadratics in exact rationals. It never calls `Surface::implicit`, a shape-side predicate or the boolean implementation to obtain a containment expectation. Primitive poses are inverted exactly. Named result boundaries use the regularised-material answers in the description, rather than applying ordinary closed-set subtraction to boundary points.

For orientation, write each primitive constraint along a normal step as `c + l t + a t²`. With `t <= 1`, choosing `t < |c| / (2 (|l| + |a| + 1))` preserves every nonzero constraint's sign. For `c = 0` and `l != 0`, `t < |l| / (2 (|a| + 1))` preserves the first-order sign. The smallest derived step is used. No fixed geometric epsilon decides a result.

`face_point` may return a contact that creates no entity. Tests check its exact surface equation and supply a rational non-contact witness when needed. The contact-only faces are entire primitive faces; quaternion rotations on a sphere/cylinder or an in-plane rational translation supply a witness. Primitive-face membership and separation from the other primitive are checked exactly. Oracle tests cover fallbacks for cases 12, 19, 20, 21, 26, 29, 34, 37 and 38. Imprint classification excludes isolated graph points.

The no-rational-point fixture has axis `(1,1,1)` and radius one. Projecting any hypothetical rational point onto `x+y+z=0` gives `(a,b,-a-b)` and `2(a²+ab+b²)=1`. Clearing denominators produces `2(A²+AB+B²)=D²`. The quadratic has even 2-adic valuation, so the left side has odd valuation and cannot be a square.

An independent Python `Fraction` script was run outside the worktree. It derives the block counts by exact cell decomposition, face-connected lumps, connected coplanar boundary patches and collinear edge merging. It confirmed all three operations in cases 1 to 9, both near-miss variants, every frame intermediate in case 42 and every union intermediate in case 44. It also checked the irrational strip width `4e-30 - 1e-60`, its listed interior witness, the circle-plane squares `16` and `27` in case 30, the rational cylinder rim, sphere perturbations, and the parity obstruction above. The curved result counts are the description's mathematical answers, not observations from an implementation.

## Speed measurements

The benchmark follows the earlier pieces' opt-in convention. Without activation it prints the command and exits before calling a piece 4 or 5 function. After both implementations exist:

```sh
CATALLAX_BOOLEAN_RUN_BENCH=1 cargo bench -p catallax-boolean --bench blind_speed
```

It prints:

- Whole-operation and individual-stage times for all pair variants, all three operations, identity placement and a fixed-seed rational placement.
- Stage counter increases, work, running maximum precision at each stage's end, marks, and fact/filter shares. Measurements use `stats()` snapshots without calling `take_stats()`. Refusals retain their code and stage.
- `classify_point` timings against result face count, with counts of Proven, ProvenCoincident, NotExact and Refused outcomes.
- Number growth after every step of cases 41 to 44 and five alternating unions/subtractions of blocks, cylinders and spheres under fixed-seed rational poses.

The blind benchmark is named `blind_speed` because the merged implementation supplies a separate `speed` benchmark. Both are retained; the implementation's benchmark is unchanged.

The growth metric is the largest observable stored rational, not the largest intermediate. Numerator and denominator maxima are separate bit counts over public rational fields, implicit coefficients, root boxes and reachable root-system coefficients. Opaque non-rational `Real` coordinates are reported separately. Hidden expression children, unrelated arena entries and temporary arithmetic allocations are not observable through this interface and are not claimed by the metric.

## Real-build verification

Merged main `f59a7a4ce16161535536723b9fb1876ca3cf17db`, including the piece 1 zero-proof speedup, into the test branch. The serial inventory ran every test against that build with a 600-second limit per test. All 77 tests compiled; the doctest command also passed. The case 19 tangent-sphere test completed in 85.691 seconds.

After applying the explicit attribute gates, the enabled suite reran successfully under the same per-test cap. Executable inventory confirmed that only the sixteen tests listed below are ignored and no tests were removed. Release Clippy passed for all targets with `-D warnings`.

| Suite | Passed | Assertion/refusal failure | Timed out |
| --- | ---: | ---: | ---: |
| Blind integration tests | 51 | 8 | 8 |
| Existing implementation unit tests | 10 | 0 | 0 |

At that review, passing blind cases ran by default and the following sixteen tests stayed ignored with their observed reasons. These are historical outcomes, not an attribution of every failure to the implementation. In particular, the normal-step and carrier-uniqueness assertions were in the blind suite's own oracle and witness helpers.

| Test | Observed reason |
| --- | --- |
| `contract::case_02_face_contact` | Inward normal-step oracle returns `None`, expected `Some(Inside)`, at `(1,0,1/2)` |
| `contract::case_03_partial_face_contact` | Outward normal-step oracle returns `None`, expected `Some(Outside)`, at `(2,1,1)` |
| `contract::case_09_near_miss_blocks` | Inward normal-step oracle returns `None`, expected `Some(Inside)`, at `(1,0,1/2)` |
| `contract::case_13_tangent_near_misses` | Exceeded the 600-second cap |
| `contract::case_14_inscribed_cylinder` | Orientation witness requires one result face on its carrier but observes two |
| `contract::case_24_steinmetz` | Unite refuses at Imprint with `boolean.imprint-ambiguous`: `no certified chart separated the boundary from the query` |
| `contract::case_29_point_tangent_cylinders` | Exceeded the 600-second cap |
| `contract::case_33_sphere_bite` | Exceeded the 600-second cap |
| `properties::cases_41_to_44_chains` | Exceeded the 600-second cap |
| `properties::cases_55_56_generated_material_identities_and_commutativity` | Exceeded the 600-second cap |
| `properties::case_57_every_scale_and_twenty_rational_placements` | Exceeded the 600-second cap |
| `properties::case_58_two_fresh_runs_and_eight_threads` | Inward normal-step oracle returns `None`, expected `Some(Inside)`, at `(1,0,1/2)` |
| `properties::case_59_low_budget_answers_are_sound_and_deadlines_are_honest` | Supported pair returns unexpected refusal `numbers.sqrt-not-separated` |
| `provenance::stage_marks_match_whole_operation_in_order_under_both_budgets` | Exceeded the 600-second cap |
| `stages::cases_51_52_imprint_classify_keep_and_direct_sew` | Inward normal-step oracle returns `None`, expected `Some(Inside)`, at `(1,0,1/2)` |
| `stages::case_53_reports_preserve_fresh_and_seeded_caller_statistics` | Exceeded the 600-second cap |

No production code or expected answer was changed to obtain these results. Timed-out aggregate tests do not prove their full placement matrices, query counts, chains, or stage checks completed. Benchmark timings are still opt-in and never gate CI.

The earlier independent exact scripts passed. Before the real build was available, a release driver also exercised the report checker with real `Decider` snapshots: fresh and seeded histories passed, while a decreasing precision peak, an incorrect counter, and cleared caller totals were rejected. That helper smoke did not execute booleans and does not substitute for the timed-out case 53 run.

## Piece 5 comparison run

The comparison at main `aea64e21dd6f6a015dd14e0bf0e55a5e53609b09` ran all 77 release tests, including every ignored body, in separate processes with the same 600-second observation cap. It recorded 62 passes, eight assertion/refusal failures and seven timeouts. Doctests passed.

`contract::case_13_tangent_near_misses` completed its whole body in 418.707 seconds and was enabled by the comparison. No expected answer or test body was changed. The comparison left the other fifteen unsuccessful bodies ignored. Their reasons distinguished test-helper defects, the Steinmetz implementation refusal, and incomplete runtime observations. A timeout is not a returned refusal or proof of a wrong answer.

The original opt-in blind benchmark had a 1,800-second observation cap. Subsequent scratch observations used 180 seconds per remaining cell and 300 seconds per growth chain. Those historical experiments are not part of this public tree.

## Post-comparison repairs

The normal-step oracle now evaluates strict regularized-result membership even when a step stays on a discarded operand boundary. An interior operand supplies a union neighbourhood; an exterior operand excludes an intersection neighbourhood; subtraction is exterior wherever its first operand is exterior or its second is interior. Two box-face interiors with opposite normals supply complementary half-neighbourhoods at a discarded union seam. Generated queries still use the original boundary-free `set_rule`, so their prescribed operand-boundary exclusions do not change.

The orientation witness no longer assumes one result face per carrier. Its replacement point must also avoid tangent contact with every other result carrier, not just the other operand's boundary. Case 59 accepts the upstream refusals `numbers.sqrt-not-separated` and, under the separately authorized low-budget scope amendment, `numbers.precision-limit`. The refinement limit actually reached supplies its own code, so a 64-bit precision limit can precede work exhaustion. None of these changes alters an assertion about a proven result.

Steinmetz failed while imprinting the second cylinder's side at query `(-3/2, 0, -1)`. The cylinder has origin `(-2, 0, 0)`, axis `(1, 0, 0)` and reference direction `(0, 0, -1)`. Its diamond phase has period four. The query is not an intersection node. The failed join is inside an ellipse trace, near phase `-2.65217391304` and axial coordinate `1.52941176471`, rather than at a singular crossing.

The old polygon builder recovered a chart translation by subtracting a Ball's rounded midpoint before and after alignment. Adding an exact period can round that midpoint, so the difference need not be a period. Applying that difference to adjacent boxes moves them into different charts by a rounding residue. Increasing precision repeats the same error. The repair carries the exact rational `4 * turns` to every box in the segment and to the closing join. Endpoint overlap, boundary separation and winding proofs remain mandatory.

Release stack samples of the timeout bodies also found algebraic elimination in `probe_path::P::excludes_zero`. A path parametrized on a surface makes that surface's defining polynomial identically zero. The adjacency proof now preserves that identity when a boundary equation is exactly the carrier polynomial, rather than rebuilding zero as a deep algebraic expression. Other equations still require the same Bernstein sign and endpoint proofs.

For near-contact path coefficients, a certified enclosure now precedes algebraic zero elimination when the fast filter cannot separate the sign. Only a ball excluding zero supplies a strict sign; unresolved coefficients retain the exact decision path. This avoids premature elimination without accepting a cached coordinate, dropping an incidence proof or changing the caller's budget.

Planar polygon winding now carries an endpoint's exact projected coordinates into the next coedge only when their vertex handles agree. It also reuses the query-coordinate differences in the crossing and between tests. The same side, crossing, segment and winding predicates still decide the result; disconnected endpoints are projected independently.

Outgoing germs on explicit lines or intersections of two planes use zero normal acceleration: their defining equations are affine, and piece 4's preceding tangent query certifies the direction. Nonaffine carriers retain the gradient, Hessian and rank proofs. This avoids reconstructing linear systems whose right-hand sides are identically zero.

Containment first tries an existing rational plane normal before the world-direction enumeration. For a capped cylinder this gives an axial ray, avoiding an unnecessary high-height quadratic on its side after tiny scaling and a large translation. This repairs the seeded `13-bite` intersection at scale `10^-28`: the captured world-X ray exhausted its root-isolation budget, while the axial side equation is a certified nonzero constant. Every intersection, forward-hit, trimmed-face, contact and parity proof remains required. Degenerate rays still try another direction; refusals still propagate unchanged.

The shared parity implementation moves an owned direction into the line and converts the submitted ray frame to exact reals only at its first hit, reusing it for later hits. Line construction and storage preserve that origin and direction. Sewing keeps its existing borrowed entry point and directed-ray semantics.

On numeric base `8fa3e04`, with the projection and affine-germ optimizations, the complete unchanged case 58 determinism body passed in 2,407.796 seconds wall time (9,767.097 seconds user CPU across its workers) and is enabled. This qualifies its two fresh runs and eight threaded runs, including every chain; it does not establish the default-suite runtime target.

On numeric base `607bbda`, the complete generated 55/56 body exposed an implementation refusal after 2,945.829 seconds, not another observation timeout. For case `29-bite` intersection, the exact point `((2-e)*3/5, (2-e)*4/5, 0)`, `e = 10^-30`, is outside the second cylinder, but `classify_point` refused with `numbers.sqrt-not-separated`. The upstream irreducible-quartic repair on `1a68cfe` makes that unchanged isolated query return `Proven(Outside)`.

Before the containment-ray change, on `1a68cfe`, the complete unchanged generated 55/56 body passed all 37,939 exact result-containment questions, topology checks and commutativity checks in 3,103.952 seconds wall time (2,946.631 seconds user CPU). Its ignore attribute was removed only after that complete pass. No generated query or expected answer was dropped; the full family remains enabled for subsequent workspace validation.

On numeric base `607bbda`, the complete unchanged stage-mark provenance body passed in 1,246.345 seconds wall time (1,233.676 seconds user CPU) and is enabled. Every default- and small-budget comparison retained mark order, whole/manual outcome agreement, counts, topology snapshots and incidence provenance.

At repair commit `13625726047923480a484f2406bbbf36193bb373`, six original bodies remained quarantined: stage cases 51/52 and 53, generated cases 55/56, placement case 57, determinism case 58, and the stage-mark provenance family. The complete stage and determinism bodies exceeded the 600-second observation cap after their oracle repair; that cap is not a mathematical failure.

The original case 57 matrix exposed a separate refusal. At scale `10^-30`, with the first pose from seed `0x570020`, case `13-bite` intersection constructs an unmarked result but refuses the transformed named inside point `(5, 5, 10 - 10^-31)`. Localization reaches `probed_region` through `geometry::on_arc`, where cyclic ordering exhausts its certified decision. The architect attributed this to piece 1 ending refinement after an unaffordable separation-bound attempt. The point must still classify `Proven(Inside)`; no piece-5 workaround or acceptance of the refusal is added.

The original literal case 57 body requested 140 pair operations and 12 chain operations for each of 61 scales and 20 poses: 185,440 booleans, before containment checks. The architect ruled that family oversized relative to the few-hundred-results and default-runtime guidance. Its approved replacement is the 456-boolean default family and an ignored 18,392-boolean qualification family. The clarified rule is: "The steps of a chain share the variant of its first operation, and `i` counts every step." Expected counts, named points and every intermediate chain check are unchanged.

At repair commit `13625726047923480a484f2406bbbf36193bb373`, rebased onto main `189fc852557a5d27ddfeb3a29c40dc3186e2952f`, `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`, and `cargo test --workspace --release -- --test-threads=1` passed. Piece 5 had 71 passing tests and six ignored bodies. All 45 contract tests passed, including Steinmetz, point-tangent cylinders with their tiny bite, and the sphere bite. Those results precede the approved case 57 restructuring and do not qualify the remaining bodies.

## Current-base shipment verification

The repair branch was rebased onto main `878f5ce010cc4187e2a04f71d6834954790be95b`, retaining its imprint decisions, trace evaluation and region-signature implementation. Main's `align` already returns the exact chart-period translation; that version is kept instead of adding a second alignment API. The containment-ray change also keeps main's independent `DIRECTION_ATTEMPTS` bound.

On repair head `ba2ac8b7cc8f02468769a0310583d588b5a26bd6`, fresh format and workspace all-targets Clippy with warnings denied passed. All fourteen enabled bodies from the comparison's fifteen unsuccessful bodies then passed complete release runs in separate processes. The generated 55/56 body checked 37,939 exact result-containment questions; determinism retained its two fresh runs and eight threaded runs. Per-body release runtimes and source fingerprints are recorded in the shipment evidence.

The actual CLI answered all three Steinmetz operations and the exact seeded `13-bite` intersection at scale `10^-28`. Each manifest was exact, unmarked and refusal-free, with the prescribed topology and surface/curve-kind counts. This verifies those four fixtures, not the deferred placement matrix or all case-13 variants.

Under the current shipment direction, the eight forge checks are the full-suite gate: format, Clippy, release archive/doc tests, four release partitions and their aggregate check. The archive inventory verifies that every enabled test belongs to exactly one partition. No second full workspace suite is run locally. The earlier local attempt on numeric base `1a68cfe` was canceled when main moved; it is not a passing full-suite run. The three recorded placement refusals and further case-13 investigation are assigned to the next fix batch.

## Focused refusal regressions

`properties.rs` runs the deferred placement cells independently of the aggregate case 57 families:

| Test | Case and operation | Scale | Common quaternion |
| --- | --- | --- | --- |
| `crossing_cylinders_tiny_rotated_union` | 23 unite | `1e-20` | `[40,-42,-1,45]` |
| `crossing_cylinders_rotated_intersection` | 23 intersect | `1/10` | `[42,66,-43,64]` |
| `crossing_cylinders_tiny_intersection` | 23 intersect | `1e-30` | identity |
| `near_tangent_cylinders_tiny_union` | 29-bite unite | `1e-30` | identity |
| `tiny_rotated_block_cylinder_cap` | 13-bite intersect | `1e-28` | `[78,44,-89,-60]` |

Each cell retains the existing exact topology counts, topology rules, material witnesses and generic rational containment query. Refusals and marks are not accepted. These focused tests do not claim completion of the aggregate placement matrix.

`coincident_sections_preserve_a_cavity_and_curve_contact` checks a three-step union, subtraction and clipping chain. Its independent geometry requires two lumps, three shells, a retained spherical cavity and a surviving cylinder/box contact. Exact queries cover the material, cavity, contact and clipping plane. This regression exercises a section that must be partitioned consistently with coincident operand edges; it does not relax the radial-order refusal.

The osculating-conic unit test independently derives third- and fourth-order terms, uses unequal tangent speeds, checks reversed comparison order and retains the duplicate-germ refusal. The actual three-dimensional contact is fourth order: on `(x-3)^2+(y-2)^2+z^2=4`, the cylinders `x^2+y^2=9` and `(x-1)^2+y^2=4` differ first in `x-3`, at coefficients `-1/96` and `-1/64` of `z^4`.

`multi_component_probe_does_not_infer_a_branch` exercises both generators of a cylinder cut by a meridian. Exact membership of the cylinder and plane must not create incidence on the other generator. Public curve-membership outcomes are checked against independent line-distance signs.
