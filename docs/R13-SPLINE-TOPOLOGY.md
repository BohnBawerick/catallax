# R13: bounded topology-search experiment

Research snapshot: 11 November 2026. Baseline:
`439aa03061fef1078bc82688463fe7f627fde4a0`.
This is partial R13 evidence, not completion of the M25 comparison or a production
intersection algorithm. No proven answer, tolerance, allowance or naming contract
changes. The M9 solver and trace-port lanes are independent.

## Sources and access limits

The methods intended by [PLAN section 12.6](PLAN.md#126-robustness-and-many-core-steps)
are identifiable from its DOI links and the [author's publication list](https://yangjieyin.github.io/publications/):

1. Jieyin Yang, Xiaohong Jia and Dong-Ming Yan, **Topology Guaranteed B-Spline
   Surface/Surface Intersection**, ACM Transactions on Graphics 42(6), 2023,
   [DOI 10.1145/3618349](https://doi.org/10.1145/3618349).
   The accessible primary material is the author's
   [SIGGRAPH Asia 2023 slides](https://yangjieyin.github.io/homepage/ACMTOG2023/SIGA23-slides.pdf),
   not the complete paper. Pages 11–17 describe Dixon implicitization, pulling
   the implicit equation onto the other parameter domain, critical-point
   elimination, tracing, inversion and clipping. The slides explicitly warn
   that an implicit surface extends beyond the original domain; a real preimage
   inside that domain is needed. Page 19 illustrates singular and tangent cases;
   examples are not general theorem hypotheses. Page 27 identifies efficiency
   as a limitation. Full degree, base-point, degeneracy and termination conditions
   could not be verified from this source.
2. Jieyin Yang and Xiaohong Jia, **A Robust and Efficient Intersection Algorithm
   for NURBS Surfaces: Handling Small Loops and Tangent Intersections**, ACM TOG
   45(5), 2026, [DOI 10.1145/3807948](https://doi.org/10.1145/3807948).
   The indexed primary abstract describes winding numbers and subdivision of a
   vector field on one parameter domain, including isolated critical points and
   nonisolated critical curves. It is insufficient to reconstruct the field or
   certify the algorithm. The [25 August 2026 erratum](https://doi.org/10.1145/3839108)
   corrects the article number from 37 to 169 according to its indexed text.
   No mathematical correction is established by the available text.

ACM returned HTTP 403 for full content. No accessible full author/preprint
version was verified; publisher restrictions were not bypassed. Consequently,
this note does **not** assert that either paper proves a theorem for arbitrary
rational patches, arbitrary degeneracies or finite-precision predicates.

The older primary reference is Patrikalakis, Maekawa and Cho,
[*Shape Interrogation for Computer Aided Design and Manufacturing*, online
section 5.8.2.3](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html).
It explains that marching needs a starting point on every branch; boundary and
collinear-normal points supply the search structure. Bounding normal pyramids
can exclude loops, whereas parallel normals cause marching difficulties.
Certified arithmetic is needed for exclusion. A loop-exclusion result alone
neither enumerates open branches nor proves disjointness.

## Comparison and certification obligations

The following obligations are our analysis of the exposed mechanisms, **not**
unverified quotations of the papers' theorems.

| Aspect | 2023 exposed mechanism | 2026 exposed mechanism | Normal-cone exclusion |
|---|---|---|---|
| Interior loops | Critical topology of pulled-back implicit equation must be complete | Winding/subdivision must locate the relevant critical structure | Can rule out loops only where certified normal sets exclude collinearity |
| Tangent/singular contacts | Multiplicities, isolated roots and branch adjacency need exact handling | Isolated versus nonisolated critical sets need distinct treatment | Failed exclusion is unresolved; parallel normals do not classify contact |
| Domain and weights | Real preimages and clipping are essential; rational denominators cannot silently vanish | Exact field and rational-denominator assumptions unavailable | Need regular charts and bounds over the whole domains, accounting for either normal orientation |
| Degeneracy | Identically zero eliminants, coincident patches and base points require explicit handling | Zero field, critical curves and cancellation require explicit handling | Zero normals or overlapping cones cannot justify exclusion |
| Termination | Root isolation and adjacency completion need hypotheses; slides insufficient | Subdivision stop/completeness conditions unavailable | A successful strict separation can stop; failed subdivision has no general completion guarantee here |
| Certified predicates | Resultants/ranks, root multiplicity/order, real inversion, domain membership, adjacency | Boundary nonvanishing, degree/winding, zero isolation and subdivision coverage, once the actual field is known | Whole-domain derivative enclosures and strict cone separation |
| Approximation | Topological equivalence does not by itself bound geometric deviation | Same distinction; no numerical deviation theorem verified | No approximation/deviation bound supplied |
| Practical cost | Exact elimination may grow degrees and coefficient sizes | Exact boundary predicates and refinement may dominate | Cheap rejection when strict bounds succeed; inconclusive regions require another method |

Neither a successful floating-point trace nor a sampled empty region supplies
these certificates. Exact rational polynomial operations are possible for the
small controls below. This establishes **no** general 20M-allowance result:
this standalone example does not measure production proof work units, run
Dixon elimination or run the 2026 field algorithm. There is no algorithm timing
ranking and no adoption recommendation yet.

## Public-case availability

The public tree at the recorded baseline contains no M25 case-set definition.
The public trace-foundation PR #388 refers to a 13-case M25 study, but does not
provide those surface inputs. Its low-level scalar/vector/box fixtures are not
that case set. Private studies, the companion repository and saved-answer
inventories were not consulted. **Mapping these results to M25 is blocked on
public case definitions.**

As a bounded fallback, the experiment reconstructs the public Rust controls in
`crates/catallax-intersect/tests/patch_fixtures/mod.rs` (`height_patch(true)`) and
plane levels in `crates/catallax-intersect/tests/patch_plane.rs`. These are not
claimed to be confirmed M25 members. No C++ candidate traces are used.

## Executable experiment and original argument

[The Rust example](../crates/catallax-features/examples/r13_topology.rs) uses
exact rationals, Cox–de Boor polynomial expansion and Bernstein sign bounds.
It implements neither published SSI algorithm. It exercises one tractable
specialization of critical-point/level reasoning and an elementary gradient
index control; the latter is **not** the unverified 2026 vector field.

The public bicubic patch has unit weights, knots
`[0,0,0,0,1,2,3,3,3,3]` in both directions, controls `(i,j,z_ij)` for
`i,j=0..5`, and `z_ij=1` exactly when both indices are in `{2,3}`.
Thus `S(u,v)=(X(u),X(v),H(u)H(v))`, where `H=N_2+N_3` and
`X=sum i*N_i`. The code verifies joins, positive `X'` on each open span,
positive interior `H`, endpoint zeros, strict increase/decrease either side
of `3/2`, and `H(3/2)=15/16`. Bernstein coefficients bound whole intervals;
these are not sample-sign assertions.

It follows that the height has a unique interior maximum `225/256` at spatial
`(5/2,5/2,225/256)`. For `0<c<225/256`, each side's inverse of `H` gives the
two halves of a single regular closed level curve, joined at their two extrema.
At the maximum the level is one point; above it the level is empty. This is a
specialized mathematical argument whose premises are checked, not a new
Catallax authoritative proof type. The boundary level zero returns `Unresolved`.

On `[1,2]^2` the height gradient points inward along all four edges, so its
nonvanishing boundary field is homotopic to the inward radial field: index +1.
This remains true for **every** plane height. Both surfaces have parallel
normals at the hill peak, so whole-patch normal exclusion cannot succeed.
Nevertheless their intersection topologies differ. Critical-point evidence
needs level/domain reasoning before it decides intersection.

Additional original controls distinguish other insufficient tests:

- Graph `F=x` against `z=0` has nonparallel constant normals but an open line:
  no closed loop does not mean no intersection.
- `F=(x^2-1)^2+y^2` on `[-1/2,3/2] × [-1,1]` has a saddle `(0,0)` and
  minimum `(1,0)`. The code verifies their stationary equations and Hessian
  signs; the explicit gradient factors exclude other roots in this rectangle.
  Their indices cancel to zero. A zero degree alone cannot prove absence.
- `F=x^2` has a critical line; isolated-zero reasoning is inapplicable.
  The linear and critical-line rows are analytic controls, not generic solvers.

## Reproduction and observations

```sh
cargo run --release -p catallax-features --example r13_topology -j2
cargo test --release -p catallax-features --example r13_topology -j2
cargo clippy -p catallax-features --example r13_topology -j2 -- -D warnings
cargo fmt --all -- --check
```

The run has nine fixed rows, no random seed, no adaptive unbounded search and
at most four concurrent copies. Each worker repeats the entire bounded control
set; equality at 1/2/4 workers tests repeatability, not parallel solver scaling.

| Control | Observed classification |
|---|---|
| Hill plane `3/10` | OneInteriorLoop |
| Tiny plane `225/256 - 1/1000000` | OneInteriorLoop |
| Tangent plane `225/256` | IsolatedContact |
| Above plane `225/256 + 1/1000000` | Disjoint |
| Boundary plane `0` | Unresolved |
| 18×18 uniform parameter samples | Maximum `1159929/1336336`, below tiny level despite its loop |
| Linear graph | No closed loop; open line |
| Cancellation graph | Index zero; two critical points; absence unresolved |
| Critical-line graph | Isolated index inapplicable |

Observed validation: the release run completed, both example tests passed,
focused Clippy passed with warnings denied, and formatting/diff checks passed.
An initial attempt to host the example in `catallax-numbers` encountered its
existing native development dependency requiring unavailable `m4`; the final
example lives in `catallax-features` and does not need that dependency. No native
oracle was run.

All hill rows report normal-rule unresolved and gradient-index +1. Two focused
unit tests check these controls and catch an endpoint-positive polynomial that
is negative inside. No production algorithms, full intersection suites or CI
runs are part of this experiment.

To finish R13: obtain the public M25 definitions and readable full primary
papers, extract their actual hypotheses, then implement and budget certified
components on those inputs. The present evidence cannot qualify either full
method or decide whether it fits Catallax's allowance.
