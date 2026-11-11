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

Keep the 2023 algebraic route as the next authoritative search prototype. Its
critical projection fits the allowance on these inputs, while rank, root
isolation, branch adjacency and clipping still need certificates. For general
spline pairs, retain a typed refusal until Catallax can certify a complete
search. Neither complete paper algorithm is qualified by this experiment.

Use the 2026 grid method only to propose critical points and starts, followed
by exact checks. Reading the full text makes this restriction firmer: Section 6
explicitly limits the guarantee to the grid resolution and acknowledges missed
loops when multiple critical points occupy one cell. Rational arithmetic does
not remove that coverage gap. The exact core fits the allowance and supplies
certified starts for the hill controls, but still refuses global coverage on
all 24 pairs. The M25 wavy pair remains unresolved. This is a proof-route
verdict, not a performance ranking of the complete algorithms.

## Sources and access

Both complete papers and their appendices were supplied and read: the
16-page 2023 paper and the 19-page corrected version of record of the 2026
paper, published on 27 August 2026. Both implemented cores are checked against
full text. The author's
[SIGGRAPH Asia slides](https://yangjieyin.github.io/homepage/ACMTOG2023/SIGA23-slides.pdf)
were also read; the supplied and downloaded copies have identical SHA-256
`c14d21da8be9fa6463d3d98e407360dfda47614b096f7d8f26a06d7c38efde4c`.

- Jieyin Yang, Xiaohong Jia and Dong-Ming Yan, *Topology Guaranteed B-Spline
  Surface/Surface Intersection*, ACM TOG 42(6), Article 211, 2023,
  [DOI 10.1145/3618349](https://doi.org/10.1145/3618349).
  Sections 4.1-4.4 specify OBB rejection, moving-plane Dixon implicitization,
  scalar critical topology, tracing and inverse-domain clipping. Section 4.3,
  Eq. 15, defines the critical resultant implemented here.
- Jieyin Yang and Xiaohong Jia, *A Robust and Efficient Intersection Algorithm
  for NURBS Surfaces: Handling Small Loops and Tangent Intersections*, ACM TOG
  45(5), Article 169, 2026,
  [DOI 10.1145/3807948](https://doi.org/10.1145/3807948).
  Equations 2-5, Section 4.1 and Definition 3 define the distance-field grid
  method; Section 4.2.2 supplies the isoparametric start strategy. Section 6
  and Appendices A-D state its coverage, regularity and resolution limits.

No paper files are committed. Raw evidence retained privately.

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

### The 2026 distance-field core

The actual field is not the height-difference gradient. With q the unique
closest point on the second surface, Eq. 2 defines
`phi = n2(q) dot (r1-r2(q))`; Eq. 3 gives
`grad(phi) = (n2(q) dot r1u, n2(q) dot r1v)`.
Appendix A distinguishes multiple closest points from a boundary closest point
whose displacement is not normal. The latter can make this signed-distance
interpretation ill-defined. Regular nonzero normals are another premise.

[The 2026 core](../crates/catallax-features/examples/r13/paper2026.rs) implements:

1. A 5 by 5 initial grid on each surface. The nearest target grid node is found
   by exact squared-distance comparisons, with deterministic tie handling.
   This implements the initial field estimate of Section 4.1.1, not a certified
   closest point on the continuous surface. Ties and selected boundary nodes
   are counted separately.
2. Eq. 3's dot products at that candidate, Eq. 5's five quadrant classes and
   Section 4.1.2's two- or four-way subdivision, with maximum depth 2. The
   unnormalized normal differs by a positive factor, so its quadrant and zero
   tests need no square root. The four-corner permutation tests of Figure 11
   record candidate indices +1 and -1. A candidate index of zero is never an
   absence certificate.
3. The four exact polynomial residuals of Eq. 4. A zero residual establishes a
   stationary pair, not that the second point is closest. No Newton iteration,
   residual tolerance or distance threshold is substituted for that distinction.
4. For the hill controls only, an independent global closest-point proof at
   the proposed peak pair, and exact start-point isolation along the two
   isoparametric curves through that peak, following Section 4.2.2.

The initial grid and depth are fixed experiment parameters, not a completeness
claim. Appendix D's general experiment uses density 15 and depth 1. This
bounded implementation does not implement general continuous closest-point
isolation, critical-curve solving, boundary intersections, tracing or a proof
that every component received a start.

#### A printed edge-rule ambiguity

Definition 3 prints `abs(chi_left-chi_right) mod 4 <= 1`, but its next sentence
says equal or adjacent quadrants. The printed formula rejects quadrants 4 and 1,
which are adjacent. This was checked on the rendered page 6, not inferred from
text extraction. The experiment runs both the literal printed test and cyclic
adjacency `min(d,4-d)<=1`, and publishes their separate invalid-edge counts.
It does not silently attribute the cyclic interpretation to the printed formula.
Neither interpretation is used as a production proof.

#### Certified peak and start witnesses

For `A.z=H(u)H(v)` and `B.z=(H(u)H(v)+c)/2`, the peak parameters are
`u=v=s=t=3/2`. All 36 controls and the four Eq. 4 residuals are checked exactly.
The normals there point along positive z. The experiment recomputes global
Bernstein enclosures `X' >= m=1` and `abs(H'') <= M=3`. With `h=15/16`,

```text
K = h*M/(4*m^2) = 45/64
B_peak.z - B(s,t).z <= K * horizontal_distance_squared
e = B_peak.z - A_peak.z
margin = 1                         when e <= 0
margin = 1 - 2*e*K                 when e > 0
```

Taylor's bound from `H'(3/2)=0` gives the height-drop bound. Hence squared
distance to any other target point exceeds the peak-pair distance by at least
`margin * horizontal_distance_squared`. A positive margin and the injective
x/y chart prove that the peak is the unique closest point. The measured margin
for the disjoint pair is `12799991/12800000`; the other three hill margins are 1.
Thus the peak field and signed gap are actual distance-field witnesses, rather
than nearest-grid estimates. Gaps are reported in the canonical frame used by
all predicates, after the checked common inverse transform.

For each loop level, the two isoparametric curves each have two starts.
The equation `h*H(t)-c=0` has one root on each strict monotone half. Twelve
exact bisections give disjoint rational isolating intervals for those roots.
The intervals define four exact root-selected parameter points; they are not
four approximate coordinates accepted at a tolerance. For the small loop the
left and right intervals are respectively
`[3069/2048,12279/8192]` and `[12297/8192,3075/2048]`.
At the contact level the single start is the exact peak, with signed gap zero.
The disjoint level has no starts on these two isoparametric curves. Its complete
empty classification still comes from the separate scalar-family proof.

A focused counterexample sets `c=100`. The peak pair still satisfies Eq. 4,
but a target boundary point is closer. The nearest-point certificate refuses
that pair. This test prevents accepting a stationary pair as a closest pair.

The wavy pair and open control fail the hill-family identity check, so no peak
or start certificate is claimed for them. On all 24 pairs the overall 2026
result is the typed refusal `r13.distance-field-coverage-unproved`, with the
measured candidate evidence and any certified witnesses retained beside it.
This replaces the former missing-text refusal with an implemented, checked
core and a specific remaining proof obligation.

#### What the full text rules out as a completeness proof

Section 5.4 and Section 6 explicitly describe lost tangent branches or small
loops when the initial grid or subdivision is too coarse. Section 6 leaves
multiple critical points in one cell for future work. Theorem 1 requires a
continuous field, a nonvanishing boundary and isolated critical points; it
cannot turn degree zero into an emptiness certificate. Appendix A's closest
point conditions and Section 6's two-dimensional overlap limitation are further
obligations. Appendix C analyzes a sufficiently small neighborhood of a critical
curve; the experiment does not assume that an arbitrary finite grid has reached
that neighborhood. Exact arithmetic alone certifies none of these coverage
premises. This is why the 2026 route is suitable for proposing candidates but
is not selected as Catallax's authoritative complete search.

### Independent height-gradient control

The earlier diagnostic remains separate. On `[1,2]^2` it computes the boundary
winding of `grad(A.z-B.z)`. Every polynomial boundary segment must stay in a
certified strict coordinate half-plane before a homotopy to its endpoint chord
is accepted. Exact oriented ray crossings give the integer winding. A boundary
zero returns `r13.boundary-field-not-separated`.

This is not the 2026 field. All four hill levels have index +1 despite their
different intersections. The original index-cancellation and critical-line
controls also remain. These show why winding evidence needs level, domain and
coverage proofs before it can classify an intersection.

## Work accounting and determinism

Each method gets a fresh, unchanged `20_000_000` operation allowance through
`with_operation_budget`. Arithmetic reserves work before execution through
`Decider::charge_work`, following `numbers::Work::exact`: rational size L is
the sum of numerator and denominator 64-bit limbs, each at least one; an add
or subtract costs `La+Lb`, and a multiply or divide costs `La*Lb`. Each
reservation adds one checkpoint. Comparisons reserve two products. This is
logical work, including no wall-clock interpretation or speed target.

Inverse transforms, basis expansion, coefficient arithmetic, Bernstein bounds,
Sylvester determinants, nearest-grid comparisons, quadrant subdivision,
stationary residuals, closest-peak certificates, start isolation, identity
checks and winding subdivision are included. Immutable fixture
construction, parsing constants, memory allocation and original nine controls
are excluded. No algorithm under comparison borrows work from another method.
An exhausted reservation aborts with `numbers.operation-budget`; its discarded
result cannot become a geometric decision. A test runs the real method with
allowance 1 to verify that behavior.

These numbers qualify only the implemented bounded tests. They provide no
20M bound for the full papers, a general SSI operation or a production Boolean.
The cone test and critical-projection costs cannot rank the complete papers.
The worker test partitions the 24 independent cases across 1, 2 and 4 threads,
then compares every ordered answer, work count and evidence record. It does not parallelize one
case or measure many-core solver scaling.

## Results

All 24 rows have identical answers, work and evidence at 1, 2 and 4 workers.
Every implemented test stays below 20M. Maximum work is 1,760,474 for normal
exclusion, 1,146,846 for the 2023 critical projection, and 1,073,095 for the
2026 core including both edge-rule interpretations and the witness checks.
A completed projection or candidate grid still leaves the M25 wavy topology
unresolved.

| Case | Normal cone | 2023 projection | Scalar family | Height-gradient winding | 2026 core |
|---|---|---|---|---|---|
| m25-wavy/base | refused; 1,760,474 | degree 15; 1,146,846 | refused; 31,494 | index 0; 100,716 | coverage refused; 1,073,095 |
| m25-wavy/moved | refused; 1,760,474 | degree 15; 1,146,846 | refused; 31,494 | index 0; 100,716 | coverage refused; 1,073,095 |
| m25-wavy/turned | refused; 1,760,474 | degree 15; 1,146,846 | refused; 31,494 | index 0; 100,716 | coverage refused; 1,073,095 |
| m25-wavy/scaled | refused; 1,760,474 | degree 15; 1,146,846 | refused; 31,494 | index 0; 100,716 | coverage refused; 1,073,095 |
| loop/base | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | coverage refused; 880,260 |
| loop/moved | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | coverage refused; 880,260 |
| loop/turned | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | coverage refused; 880,260 |
| loop/scaled | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | coverage refused; 880,260 |
| contact/base | refused; 1,708,176 | degree 12; 782,472 | IsolatedContact; 31,458 | index 1; 100,608 | coverage refused; 875,928 |
| contact/moved | refused; 1,708,176 | degree 12; 782,472 | IsolatedContact; 31,458 | index 1; 100,608 | coverage refused; 875,928 |
| contact/turned | refused; 1,708,176 | degree 12; 782,472 | IsolatedContact; 31,458 | index 1; 100,608 | coverage refused; 875,928 |
| contact/scaled | refused; 1,708,176 | degree 12; 782,472 | IsolatedContact; 31,458 | index 1; 100,608 | coverage refused; 875,928 |
| small-loop/base | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | coverage refused; 880,824 |
| small-loop/moved | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | coverage refused; 880,824 |
| small-loop/turned | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | coverage refused; 880,824 |
| small-loop/scaled | refused; 1,708,176 | degree 12; 782,472 | OneInteriorLoop; 31,458 | index 1; 100,608 | coverage refused; 880,824 |
| disjoint/base | refused; 1,708,176 | degree 12; 782,472 | Disjoint; 31,458 | index 1; 100,608 | coverage refused; 876,378 |
| disjoint/moved | refused; 1,708,176 | degree 12; 782,472 | Disjoint; 31,458 | index 1; 100,608 | coverage refused; 876,378 |
| disjoint/turned | refused; 1,708,176 | degree 12; 782,472 | Disjoint; 31,458 | index 1; 100,608 | coverage refused; 876,378 |
| disjoint/scaled | refused; 1,708,176 | degree 12; 782,472 | Disjoint; 31,458 | index 1; 100,608 | coverage refused; 876,378 |
| open-control/base | no closed loop; 1,683,876 | zero derivative; 653,460 | OneOpenArc; 32,550 | index 0; 100,632 | coverage refused; 829,236 |
| open-control/moved | no closed loop; 1,683,876 | zero derivative; 653,460 | OneOpenArc; 32,550 | index 0; 100,632 | coverage refused; 829,236 |
| open-control/turned | no closed loop; 1,683,876 | zero derivative; 653,460 | OneOpenArc; 32,550 | index 0; 100,632 | coverage refused; 829,236 |
| open-control/scaled | no closed loop; 1,683,876 | zero derivative; 653,460 | OneOpenArc; 32,550 | index 0; 100,632 | coverage refused; 829,236 |

Normal refusals use `r13.normal-cones-not-separated`. The scalar refusal is
`r13.scalar-topology-not-implemented`. Projection refusals use
`r13.zero-critical-derivative`. The 2026 code is
`r13.distance-field-coverage-unproved`. Neither a projection degree nor a
candidate index is an intersection verdict.

The following 2026 evidence is identical for the four copies of each family.
All discrete nearest-node tie counts are zero. This does not prove continuous
nearest-point uniqueness. Invalid edges are counted per leaf cell, so shared
edges may appear twice.

| Family | Source samples | Boundary target nodes selected | Exact grid Eq. 4 zeros | Printed / cyclic invalid edges | Cyclic candidate +1 / -1 cells | Certified starts |
|---|---|---|---|---|---|---|
| m25-wavy | 50 | 24 | 0 | 15 / 1 | 1 / 1 | unproved |
| loop | 31 | 18 | 17 | 2 / 0 | 0 / 0 | 4 |
| contact | 31 | 18 | 17 | 2 / 0 | 0 / 0 | 1 |
| small-loop | 31 | 18 | 17 | 2 / 0 | 0 / 0 | 4 |
| disjoint | 31 | 18 | 17 | 2 / 0 | 0 / 0 | 0 |
| open-control | 25 | 10 | 0 | 0 / 0 | 0 / 0 | unproved |

The hill rows include stationary boundary pairs as well as the peak. Those
counts do not enumerate isolated interior critical points. The wavy pair has
one positive and one negative candidate-index cell, but no exact grid Eq. 4
zero and no certified starts. No claim of two actual critical points follows
from those sampled indices.

Validation: the optimized example completed and ten focused tests passed.
They cover all 24 cases at 1/2/4 workers, every wavy resultant against the
kernel reference, the edge-rule discrepancy, stationary-versus-closest
counterexample, boundary zeros and budget exhaustion. Local runs use
`rustc -O` against the Cargo-built release numeric library. Focused Clippy
and packaged release Clippy pass with warnings denied. Formatting and the
frozen-input generator comparison are recorded with the private evidence.
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
