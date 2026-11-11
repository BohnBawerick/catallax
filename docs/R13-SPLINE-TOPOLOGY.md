# R13: exact spline topology-search experiment

Research date: 11 November 2026. This extends the original R13 note and Rust
controls from `dot/r13-spline-topology-comparison`, merged before this work.
Only the opt-in example and this note change. No production intersection,
answer, tolerance, work allowance or replay fixture changes.

## Verdict for the spline-with-spline stage

Use exact normal-cone exclusion as a first test. It is a sufficient test for
absence of closed loops, and must never stand in for a complete intersection
search. Follow it with exact scalar topology where a shared injective graph
chart can be proved. This experiment decides the loop, small-loop, isolated
contact and disjoint controls that way. It refuses the M25 wavy pair.

Do not select either complete paper algorithm for production from this evidence.
For general spline pairs, retain a typed refusal until Catallax can certify a
complete search. The 2023 algebraic route is the next bounded prototype to
investigate. Its critical projection fits the allowance on these inputs, and
its full text makes the remaining rank, root-separation and clipping gaps
concrete. The 2026 winding route remains a research candidate.
Its exact field and treatment of critical curves need full-text review before
implementation. This is a decision about the available proof route, not a
performance ranking between unimplemented algorithms.

## Sources and access

The complete 16-page 2023 paper, including its appendices, was supplied and
read before this note was finished. Its Section 4.3, Eq. 15, confirms the
critical-point resultant implemented here. The author's
[SIGGRAPH Asia slides](https://yangjieyin.github.io/homepage/ACMTOG2023/SIGA23-slides.pdf)
were also read; the supplied and downloaded copies have identical SHA-256
`c14d21da8be9fa6463d3d98e407360dfda47614b096f7d8f26a06d7c38efde4c`.

- Jieyin Yang, Xiaohong Jia and Dong-Ming Yan, *Topology Guaranteed B-Spline
  Surface/Surface Intersection*, ACM TOG 42(6), Article 211, 2023,
  [DOI 10.1145/3618349](https://doi.org/10.1145/3618349).
  Sections 4.1-4.4 specify OBB rejection, moving-plane Dixon implicitization,
  scalar critical topology, tracing and inverse-domain clipping.
- Jieyin Yang and Xiaohong Jia, *A Robust and Efficient Intersection Algorithm
  for NURBS Surfaces: Handling Small Loops and Tangent Intersections*, ACM TOG
  45(5), 2026, [DOI 10.1145/3807948](https://doi.org/10.1145/3807948).
  This method is **not checked against full text**. The
  [publisher-deposited abstract](https://api.crossref.org/works/10.1145/3807948)
  describes winding numbers and subdivision of a field on one parameter domain
  to find starting points, including small loops and tangent branches. It does
  not specify enough of the field to reconstruct that algorithm. An accompanying
  supplied method summary was read; information from it is from summary, not
  full text, and supplies no additional certified predicate here.

Initial publisher requests returned challenge pages, and the requested browser
attempt failed with `BRIDGE_NOT_READY`. The 2023 full text subsequently arrived
through the supplied reading material. No paper file is committed. Raw evidence
retained privately.

### What the 2023 full text changes

The paper is not a purely floating-point algorithm. Section 4.3 explicitly uses
exact computation to isolate real roots of its critical resultant. It also
contains practical steps that need stronger certificates in Catallax:

- Proposition 1 distinguishes full-rank Dixon determinants without base points
  from implicit factors of maximal nonzero minors with base points. Remark 1
  estimates polynomial rank by generic evaluations. A finite set of successful
  evaluations is not an upper-rank proof. Exact rank, factor and extraneous-root
  checks remain necessary before adopting that implicit equation.
- Theorem 4.1 and Remark 2 require real parameter preimages in the represented
  domain, with a nonzero homogeneous point. A point on the implicit extension
  is not automatically on the patch. Remark 3 uses selected matrix rows followed
  by a distance threshold to filter candidate inverses; that threshold cannot
  replace an exact membership proof.
- For rational critical abscissae, Section 4.3 uses a perturbation of `1e-6` in
  its experiments. Catallax needs a certified root-separation choice instead.
  Section 5.5 separately handles vertical components, isolated points and whole
  boundary curves. A zero resultant or missing neighboring roots is not empty.
- Section 4.4 samples boundary curves with step `1/N` to find splitting points
  before representative-point clipping. Catallax needs complete boundary
  intersections and certified branch membership to make this a proof.
- Table 1 and Appendix B use geometric tolerances and numerical tracing. The
  topology analysis does not itself certify whole-curve geometric deviation.

Appendix A exhibits a pulled-back curve of bidegree `(18,24)`. Our shared-chart
bicubics avoid that degree growth and all Dixon construction. Their measured
resultant work therefore does not bound the general algorithm. These are our
certification obligations after reading the paper, rather than a claim that
its practical implementation already provides Catallax's proof contract.

The older rule comes from Patrikalakis, Maekawa and Cho,
[*Shape Interrogation*, section 5.8.2.3](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node109.html).
It describes bounding normal pyramids for loop exclusion and the need for exact
or interval arithmetic. Marching still needs a starting point on every branch.
A failed cone test proves nothing about presence, absence or type of a meeting.

## Public inputs and their relation to M25

All cases have unit weights, degree `(3,3)`, the knot vector
`[0,0,0,0,1,2,3,3,3,3]` in each parameter, and a 6 by 6 control net.
Parameters range over `[0,3]^2`. There are nine polynomial knot cells per patch.

M25's scope study includes plane, tilt, rational, tight, loop, tangent,
neartangent, tangentline and cylcut, and one spline-with-spline pair described
as two wavy surfaces. Its appendix supplies the constructor. The present
`m25-wavy` fixture reproduces that constructor, rather than substituting the
hill controls for it. Control points have `x=i`, `y=j` and

```text
z = amp * sin(0.9*i + shift) * cos(0.7*j - shift) + 0.05*i
first:  amp = 1.0,  shift = 0.3
second: amp = -0.8, shift = 1.1
```

The [frozen fixture](../crates/catallax-features/examples/r13/wavy_fixture.rs)
contains the binary64 bits from evaluating that C++ expression with
`-O0 -ffp-contract=off`. The public
[generator](../crates/catallax-features/examples/r13/freeze_wavy.cpp) reproduces
the formatted fixture byte for byte on the checked compiler and math library. The experiment converts those bits
to exact rationals; it does not call trigonometric functions at runtime. This
makes its input portable. The original study's serialized control bytes were
unavailable, so byte identity with that earlier run is not claimed. Its reported
SMLib result is two curves, not a Catallax proof or a rerun of the native oracle.

For the other required families, let `H=N2+N3`, `X=sum(i*Ni)` and
`A=(X(u),X(v),H(u)H(v))`. Let
`B=(X(u),X(v),(H(u)H(v)+c)/2)`. Both patches are curved splines.
Their intersection is the level `H(u)H(v)=c`.

| Family | Exact c | Known topology |
|---|---|---|
| loop | `3/10` | One interior loop |
| contact | `225/256` | One isolated contact |
| small-loop | `225/256 - 1/1000000` | One small interior loop |
| disjoint | `225/256 + 1/1000000` | Empty |

An additional open control has `B.z=A.z+10*(X(u)-5/2)`. Its intersection is one
open arc. This makes a successful loop-exclusion result observable without
claiming that it proves disjointness.

Each family has a base copy, translation `(17/3,-11,2/7)`, rotation
`[[3/5,0,4/5],[0,1,0],[-4/5,0,3/5]]`, and uniform scale `7/3`.
There are 24 spline pairs in total, including four M25 wavy copies.
The experiment physically transforms both control nets, verifies orthogonality
and positive scale, then applies the rational inverse before each predicate.
This canonicalization is charged. These copies test common rigid motion and
scale invariance; they do not test independently moving one surface.

## Implemented tests and their limits

[The opt-in example](../crates/catallax-features/examples/r13_topology.rs) keeps
the original nine controls and adds
[the spline experiment](../crates/catallax-features/examples/r13/spline_search.rs).
No research result becomes a production certificate.

### Exact normal-cone exclusion

The example expands Cox-de Boor bases on all knot cells, differentiates the
spatial patches, forms `N=Su cross Sv`, and converts each component to tensor
Bernstein coefficients. Their convex hulls enclose the whole normal ranges.
It considers every pair of the nine cells on each surface, including cells
that do not meet. It accepts only if one fixed component of `Na cross Nb` has
one strict sign over all 81 normal-box pairs.

This is a conservative implementation of the older exclusion rule. The strict
cross-component test excludes both parallel and antiparallel normals without
normalizing them or taking square roots. Along any regular intersection its
oriented tangent then has a strictly monotone spatial coordinate, so a closed
component is impossible. A zero normal or overlapping bounds cannot pass.
Loose boxes may refuse even when sharper cones could separate.

A successful answer is `NoClosedLoop`, not a classification of all branches.
Boundary enumeration, component connection and geometric error bounds remain
separate obligations. The required hill pairs have parallel normals at the
peak, even when they are disjoint, so refusal there is expected.

### The 2023 critical projection and a scalar-family control

The shared x and y control identities and `X'>0` imply that equal spatial
points have the same parameters. For this case set, the scalar equation is
therefore exactly `g(u,v)=A.z-B.z`, with no implicitization or inverse clipping.
The example constructs, on all nine knot cells,
`R(u)=Res_v(g,partial_v g)`, as in the slides' page 15 and the paper's Eq. 15.
It computes a Sylvester determinant over rational polynomials in u. The matrix
is at most 5 by 5. A bounded Laplace expansion keeps each coefficient operation
interruptible under the existing allowance. The kernel's public resultant is
unbounded, so it is used only as an independent test reference. All nine wavy
resultants agree exactly with that reference; a separate quadratic identity
checks the determinant's sign and coefficients.

A successful `CriticalProjection` certifies only construction of these nonzero
projection polynomials. It does not isolate their roots, remove projection
artifacts, connect branches, or classify intersections. An identically zero
critical derivative returns `r13.zero-critical-derivative`; a zero projection
polynomial returns `r13.nonisolated-critical-projection`. In particular, the
open control refuses this fixed projection direction even though its topology
is independently known. This experiment does not implement Section 5.5's
special treatment or a change of projection direction.

A separate scalar-family control checks the height-difference identities over
all 36 controls. It proves endpoint values, joins, positivity and strict
increase then decrease of H by exact Bernstein bounds. Thus H has its only
interior maximum `15/16` at `3/2`, and the height product has maximum `225/256`.
Its two inverse branches establish a single closed level for each interior
positive level below that maximum. The maximum itself is an isolated point;
levels above it or below zero are empty. The boundary level refuses.

The open control's verified difference is linear in X, so strict monotonicity
of X gives one open arc. Other scalar equations, including M25's wavy pair,
return `r13.scalar-topology-not-implemented`. The scalar-family classification
is an original specialized proof, not the paper's general topology algorithm.
The 2023 resultant component is checked against full text; the complete paper
algorithm is not implemented or qualified.

### Gradient winding control, separate from the 2026 method

On the central parameter cell `[1,2]^2`, the example computes the boundary
winding of `grad(A.z-B.z)`. Each polynomial boundary segment must lie in a
certified strict coordinate half-plane. Otherwise it subdivides, up to depth
12 or the operation allowance. This certifies a homotopy to the endpoint
polygon that avoids zero. Exact oriented ray crossings then give its integer
winding. A boundary zero returns `r13.boundary-field-not-separated`.

This field is chosen for the experiment. It is not asserted to be the 2026
paper's field. Its index is evidence about critical points inside that one
cell, not about intersection topology on the complete domain. In particular,
all four hill levels have index +1 despite their different intersections.
The original index-cancellation and critical-line controls remain to show
why degree zero and isolated-zero logic are insufficient.

The 2026 paper method itself returns `r13.paper-field-unverified` with zero
algorithm work. The complete 2026 algorithm is not checked against full text.
Zero here means unimplemented, not fast or free.

## Work accounting and determinism

Each method gets a fresh, unchanged `20_000_000` operation allowance through
`with_operation_budget`. Arithmetic reserves work before execution through
`Decider::charge_work`, following `numbers::Work::exact`: rational size L is
the sum of numerator and denominator 64-bit limbs, each at least one; an add
or subtract costs `La+Lb`, and a multiply or divide costs `La*Lb`. Each
reservation adds one checkpoint. Comparisons reserve two products. This is
logical work, including no wall-clock interpretation or speed target.

Inverse transforms, basis expansion, coefficient arithmetic, Bernstein bounds,
Sylvester determinants, identity checks and winding subdivision are included. Immutable fixture
construction, parsing constants, memory allocation and original nine controls
are excluded. No algorithm under comparison borrows work from another method.
An exhausted reservation aborts with `numbers.operation-budget`; its discarded
result cannot become a geometric decision. A test runs the real method with
allowance 1 to verify that behavior.

These numbers qualify only the implemented bounded tests. They provide no
20M bound for the full papers, a general SSI operation or a production Boolean.
The cone test and critical-projection costs cannot rank the complete papers.
The worker test partitions the 24 independent cases across 1, 2 and 4 threads,
then compares every ordered answer and work count. It does not parallelize one
case or measure many-core solver scaling.

## Results

All 24 rows have identical answers and work at 1, 2 and 4 workers. Work is
shown after each result. Every implemented test stays below 20M. The largest
normal test costs 1,760,474 units; the largest critical projection costs
1,146,846. A completed projection still leaves the M25 wavy topology unresolved.

| Case | Normal cone | 2023 projection | Scalar family | Gradient winding | 2026 method |
|---|---|---|---|---|---|
| m25-wavy/base | refused; 1,760,474 | degree 15; 1,146,846 | refused; 31,494 | index 0; 100,716 | unverified; 0 |
| m25-wavy/moved | refused; 1,760,474 | degree 15; 1,146,846 | refused; 31,494 | index 0; 100,716 | unverified; 0 |
| m25-wavy/turned | refused; 1,760,474 | degree 15; 1,146,846 | refused; 31,494 | index 0; 100,716 | unverified; 0 |
| m25-wavy/scaled | refused; 1,760,474 | degree 15; 1,146,846 | refused; 31,494 | index 0; 100,716 | unverified; 0 |
| loop/base | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | unverified; 0 |
| loop/moved | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | unverified; 0 |
| loop/turned | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | unverified; 0 |
| loop/scaled | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | unverified; 0 |
| contact/base | refused; 1,708,176 | degree 12; 782,472 | IsolatedContact; 31,458 | index 1; 100,608 | unverified; 0 |
| contact/moved | refused; 1,708,176 | degree 12; 782,472 | IsolatedContact; 31,458 | index 1; 100,608 | unverified; 0 |
| contact/turned | refused; 1,708,176 | degree 12; 782,472 | IsolatedContact; 31,458 | index 1; 100,608 | unverified; 0 |
| contact/scaled | refused; 1,708,176 | degree 12; 782,472 | IsolatedContact; 31,458 | index 1; 100,608 | unverified; 0 |
| small-loop/base | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | unverified; 0 |
| small-loop/moved | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | unverified; 0 |
| small-loop/turned | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | unverified; 0 |
| small-loop/scaled | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | unverified; 0 |
| disjoint/base | refused; 1,708,176 | degree 12; 782,472 | Disjoint; 31,458 | index 1; 100,608 | unverified; 0 |
| disjoint/moved | refused; 1,708,176 | degree 12; 782,472 | Disjoint; 31,458 | index 1; 100,608 | unverified; 0 |
| disjoint/turned | refused; 1,708,176 | degree 12; 782,472 | Disjoint; 31,458 | index 1; 100,608 | unverified; 0 |
| disjoint/scaled | refused; 1,708,176 | degree 12; 782,472 | Disjoint; 31,458 | index 1; 100,608 | unverified; 0 |
| open-control/base | no closed loop; 1,683,876 | zero derivative; 653,460 | OneOpenArc; 32,550 | index 0; 100,632 | unverified; 0 |
| open-control/moved | no closed loop; 1,683,876 | zero derivative; 653,460 | OneOpenArc; 32,550 | index 0; 100,632 | unverified; 0 |
| open-control/turned | no closed loop; 1,683,876 | zero derivative; 653,460 | OneOpenArc; 32,550 | index 0; 100,632 | unverified; 0 |
| open-control/scaled | no closed loop; 1,683,876 | zero derivative; 653,460 | OneOpenArc; 32,550 | index 0; 100,632 | unverified; 0 |

Normal refusals use `r13.normal-cones-not-separated`. The scalar refusal is
`r13.scalar-topology-not-implemented`. Projection refusals in this set use
`r13.zero-critical-derivative`. The 2026 method returns
`r13.paper-field-unverified`; its zero work records an unimplemented method.
Neither a projection degree nor a winding index is an intersection verdict.

The original nine controls remain: they distinguish isolated contact, an
interior loop, disjointness, an open line, index cancellation, a critical
line, and the 324-sample grid that misses the small hill loop.

Validation: the optimized executable completed; seven focused tests passed,
including all nine wavy resultants against the kernel reference. These local
runs used `rustc -O` against the built release `catallax-numbers` library.
Focused `clippy-driver` and packaged `cargo clippy --release` both passed
with warnings denied. Workspace formatting and the fixture-generator byte
comparison passed. The packaged Cargo test build was
stopped while compiling unchanged dependencies; it is not claimed as passing.

## Reproduction

```sh
cargo run --release -p catallax-features --example r13_topology -j2
cargo test --release -p catallax-features --example r13_topology -j2
cargo clippy --release -p catallax-features --example r13_topology -j2 -- -D warnings
cargo fmt --all -- --check
```

No native oracle, full paper algorithm, production intersection suite or replay
corpus is claimed by these commands. The 477 saved replay answers are untouched;
this example does not rewrite them. Raw evidence retained privately.
