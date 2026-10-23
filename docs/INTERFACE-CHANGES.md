# Interface changes

An append-only log of every change to a public interface: a type, field, function signature, trait, outcome or refusal code, budget default, or file format that crosses a piece boundary or the judge boundary.

Rules:

- Public interface changes require review and a migration for affected callers.
- Every change gets an entry here in the same pull request that makes it, before or with the code.
- Entries are never edited or removed. A correction is a new entry that names the one it corrects.
- Name the affected pieces in the pull request so their callers can migrate together.

Entry format:

```
## YYYY-MM-DD: short title

Pieces: 1, 2, ... (or "judge interface", "all")
Change: what changed, in one or two sentences.
Why: the reason.
Migration: what a worker on an affected piece must do.
```

## 2026-10-18: initial interfaces

Pieces: all, judge interface
Change: first version of `docs/DESIGN.md`, `docs/ARCHITECTURE.md`, the seven piece descriptions, `docs/JUDGE-INTERFACE.md` with question format `exact-brep-question/1`, manifest format `exact-brep-manifest/1` and the `exact-brep` program, and the initial crate interfaces. The crate names in this public tree use the `catallax-` prefix. Budget default: `max_precision_bits = 16384`, `max_work = 100_000_000`.
Why: kickoff of the skeleton.
Migration: none; this is the starting point.

## 2026-10-18: question and manifest formats move to version 2

Pieces: judge interface, cli
Change: the question format is now `exact-brep-question/2` and the manifest format `exact-brep-manifest/2`, carrying the ND1 to ND11 decisions below. `catallax_cli::QUESTION_FORMAT` and `catallax_cli::MANIFEST_FORMAT` change to match. Version 1 is withdrawn.
Why: the decisions change the default deadline, when a STEP file is written and how `result` counts touching lumps, and `JUDGE-INTERFACE.md` requires a version bump for any format change.
Migration: the `catallax-cli` creator writes and reads only version 2; `about` lists only `exact-brep-question/2`.

## 2026-10-18: ND1 rotations stay rational quaternions

Pieces: judge interface, cli
Change: a pose rotation is a quaternion `[w, x, y, z]` of exact numbers; integer quaternions are the usual form. Rotation matrices are not accepted.
Why: every rotation matrix with rational entries comes from an integer quaternion, so one form reaches every exact pose and the kernel never has to check a matrix for orthogonality.
Migration: none.

## 2026-10-18: ND2 numbers are decimal or p/q strings

Pieces: judge interface, cli
Change: question numbers use piece 1's grammar, a decimal with an optional exponent or a fraction `p/q`, and nothing else. Every number in the manifest is a decimal in the notation of `to_decimal_bounds`, never a fraction.
Why: both forms parse exactly; writing only decimals keeps the manifest readable by any decimal parser.
Migration: `render_manifest` writes bounds, mark bounds and `max_deviation` with `to_decimal_bounds`.

## 2026-10-18: ND3 an empty answer is answered, with no STEP file

Pieces: 7, judge interface, cli
Change: an empty answer has status `answered`, `result.kind` `empty` with zero counts, volume and area exactly 0, centroid `null`, `step` `null` and no `answer.step`. A stale STEP file is removed.
Why: an empty answer can be reached through not-exact decisions (a near miss treated as touching), so it needs the same `exact` flag and marks as any answer, which a separate status would hide.
Migration: the `catallax-cli` creator skips `write_step` for an empty result. Piece 7's `write_step` still writes a valid file for an empty body when called.

## 2026-10-18: ND4 touching material is separate lumps sharing no entity

Pieces: 3, 5, 7, judge interface
Change: where a result's material touches itself, the two sides share no face, edge or vertex. Pieces that touch only along curves or at points are separate lumps, each its own `MANIFOLD_SOLID_BREP`; a connected piece touching itself at isolated points is one lump; a connected piece touching itself along a curve is refused with `boolean.self-touching-shell`. This replaces the rule that touching lumps may share edges and vertices.
Why: lumps with their own entities export as self-contained manifold solids, and `result` counts are then the sum of the lumps' counts; a non-manifold solid has no `MANIFOLD_SOLID_BREP` form.
Migration: piece 3 and piece 5 creators never share topology between lumps; piece 7 writes each lump's entities on their own.

## 2026-10-18: ND5 export deviation, B-spline branches, no pcurves

Pieces: 7, judge interface
Change: `max_deviation` is at most `curve_tolerance`, default 1e-7 mm, unless the answer carries `step.tolerance-not-reached`, which piece 7 also adds when rounding at `significant_digits` alone exceeds the tolerance. Branches with no closed form are one non-rational `B_SPLINE_CURVE_WITH_KNOTS` each. Edges carry no pcurves.
Why: a pcurve is a second approximation per face that would need its own proven distance, and readers build it from the 3D curve and the analytic surface.
Migration: piece 7 writes no `PCURVE`, `SURFACE_CURVE` or `SEAM_CURVE`, and checks rounding against the tolerance.

## 2026-10-18: ND6 maximal faces and edges

Pieces: 5, judge interface
Change: the kernel promises maximal faces and edges: no edge between two faces on the same surface with the same sense, no vertex joining exactly two edges on the same curve that bound the same two faces, and a ring edge where an edge covers a whole closed curve alone. This entry's original marked-merge rule is superseded by the later `boolean.merge-ambiguous` rule.
Why: the manifest counts then match the counts maths gives, and the pieces merged are proven or marked like every other decision.
Migration: piece 5's sew removes redundant vertices as well as redundant edges.

## 2026-10-18: ND7 analytic carriers stay analytic in STEP

Pieces: 7, judge interface
Change: every face surface is written as its analytic STEP entity and the file never holds a B-spline surface; lines, circles and ellipses stay analytic.
Why: this was already piece 7's rule; it is now stated in the judge contract.
Migration: none.

## 2026-10-18: ND8 questions carry no labels

Pieces: judge interface, cli
Change: the question `id` and all names are opaque; questions carry no pose, family or cell label and the kernel reads none. The examples now use opaque ids. This confirms the judge's provisional behaviour.
Why: an answer must come from the geometry alone.
Migration: none.

## 2026-10-18: ND9 running environment, limits and exit codes

Pieces: judge interface, cli
Change: the default `deadline_ms` becomes 50,000, counted from the program's start. The judge's outer limits are 60 s and 4 GiB; the kernel does not measure its memory. The program needs an empty sandbox only. Reruns of one build are byte-identical except `stats` and deadline refusals. A refusal exits 0; every nonzero exit, signal, timeout or missing manifest is a crash.
Why: a 50 s deadline leaves time to write the refusal manifest inside the judge's 60 s limit, which a 60 s deadline would race.
Migration: the `catallax-cli` creator applies the new default and sets the decider's deadline from the program's start.

## 2026-10-18: ND10 properties enclose the built body, never the STEP copy

Pieces: 6, judge interface
Change: `properties` intervals enclose the exact properties of the built body. They contain the true answer's values unless a mark from outside pieces 6 and 7 is present; marks with codes `measure.*` or `step.*` never change the body. This accepts the judge's provisional behaviour for answers without such marks.
Why: a mark from a not-exact decision means the built body may differ from the true answer, so its intervals cannot promise the true values.
Migration: none for piece 6; its enclosures were already of the built body.

## 2026-10-18: ND11 seams and poles in the operand's frame

Pieces: 2, 3, 7, judge interface
Change: every periodic surface of an operand has its seam on the local half-plane `y = 0, x > 0`, a sphere's poles lie on the local `z` axis, and a torus's tube seam is its outer equator. Booleans keep operand surfaces, so answers keep those seams, and STEP placements show them in the same place. This confirms the judge's provisional behaviour.
Why: this was already the builders' and piece 2's parametrisation; it is now stated in the judge contract so "along the seam" poses can aim at it.
Migration: none; piece 3's builders already use `ref_dir = +x` and `axis = +z` in the local frame.

2026-10-18: Piece 1, case 7's remainder rule corrected to floor semantics by architect decision; `q = floor(a / d)` and a non-zero remainder has the divisor's sign. No signature change.

## 2026-10-19: increase the piece 1 default work budget

Pieces: 1
Change: `BudgetPolicy::default().max_work` increases from `100_000_000` to `1_000_000_000`. `max_precision_bits` remains `16384`.
Why: the architect approved the retune after filtered Descartes isolation exhausted the old limit on degree-64 Chebyshev polynomials. Speed is measured, not gated; the work limit still bounds each decision.
Migration: callers using the default receive the larger work budget. Explicit policies are unchanged.

## 2026-10-19: piece 2 adds surface derivatives and a curve's parameter origin

Pieces: 2, 4, 6, 7
Change: `Surface::derivatives(&self, u: &Ball, v: &Ball) -> [[Ball; 3]; 2]` encloses `[dS/du, dS/dv]`, and `Curve::parameter_origin(&self, store: &GeomStore) -> Point3` gives the exact point at `t = 0` (a circle's or ellipse's point in its reference direction, an intersection curve's witness).
Why: piece 6 integrates over surface parametrisations and needs their derivatives; pieces 4 and 7 start closed curves and place ring-edge vertices at the parameter origin, which needs one exact definition.
Migration: piece 2 creator implements both. Pieces 4, 6 and 7 use them instead of deriving their own.

## 2026-10-19: piece 2 fixes its exact forms and validation order

Pieces: 2, 3, 4
Change: `Surface::implicit` and `Curve::equations` return exactly the polynomials of the piece 2 tables, unscaled and in the stated order; a line's equations are the non-zero components of `(x - origin) x dir` (two or three) instead of "two planes". An ellipse now needs `0 < minor_sq < major_sq`; equal semi-axes are refused with `shapes.ellipse-axes` and must be a `Circle`. Constructors check zero vectors, then scalar ranges, then perpendicularity, and return the first failure. `parameters_of` uses the stated formulas, with `u` in `[0, 2 pi)` so its jump is on the seam. `Rotation` equality compares the matrix.
Why: blind tests need one exact answer per call, and one representation per circle keeps `same_surface` and piece 4's closed forms unambiguous.
Migration: piece 4 builds an ellipse only when its axes differ and a circle otherwise. Piece 3 builders validate before calling piece 2 constructors, so they meet none of these refusals.

## 2026-10-19: piece 2 decisions state their outcomes when they cannot prove

Pieces: 2, 4, 5
Change: a `NotExact` from `side`, `side_of`, `point_on_curve` or `same_point` carries a proven distance bound in model units (callers mark it with `BoundKind::Distance`); on an intersection curve it bounds the distance to each of the two surfaces, not to the curve. `point_on_curve` and `same_point` decide equations or coordinates one at a time in a fixed order and stop at the first proven non-zero or refused one. One from `normals_agree` carries the dot-product bound (`BoundKind::Value`). `same_surface` uses exact rational arithmetic without the `Decider` and is never `NotExact` or `Refused`. A fact-settled decision calls only `note_fact_decision` and answers after the deadline. `side` answers about the implicit's zero set, so a cone's opposite nappe is coincident. The decisions take `&self` and never record facts.
Why: piece 1's comparison showed that blind tests need the exact outcome of every unproven case, and marks from these decisions must have a geometric meaning.
Migration: callers of these decisions build marks with the stated `BoundKind`. A caller that needs the cone's nappe tests `(p - apex) . axis` itself. A caller that wants a proven zero kept records it with origin `Proof`.

## 2026-10-19: piece 3 fixes its builders' geometry and refusal order

Pieces: 3, 5, 6, 7
Change: the piece 3 description now gives each builder's exact surfaces, curves, points, face order, loop cycles and coedge senses in the local frame, the cone surface's apex, axis sign and slope for both frustum directions, the facts each builder must record (with origin names `topology.block` and so on), and how many entities each adds. Builders check `topology.non-positive-size` before `topology.cone-is-cylinder` and `topology.torus-self-intersecting`, and a refused builder adds nothing to the model, the store or the fact book. No signature changes.
Why: blind tests compare builder output field by field, and later pieces read the builders' seams, senses and facts.
Migration: none for signatures. Piece 5 can rely on the listed facts for every operand edge and vertex.

## 2026-10-19: piece 3 check gains fault codes and entity-sharing rules

Pieces: 3, 5
Change: `Model::check` reports eleven fault codes: the earlier `topology.loop-not-closed`, `topology.edge-use`, `topology.euler`, `topology.disconnected-shell` and `topology.missing-incidence`, plus new `topology.body-kind`, `topology.empty-part`, `topology.shared-entity`, `topology.edge-ends`, `topology.unbounded-face` and `topology.vertex-loop`. Missing incidence now also covers an edge's end points on its curve and a vertex loop's point on its surface. The Euler count is stated exactly. Bodies may share entities with other bodies; within one body no face, edge or vertex is used by two shells. `owner_of_coedge` gives the first face added that holds the coedge.
Why: the earlier list could not report a vertex loop alone in a face, a loopless plane face or an entity shared between shells, and ND4 needs the last one checked.
Migration: piece 5 records `Incident` facts for every vertex it adds on an edge's curve (piece 2 does this for `Meet` points), so its results pass `check`. Piece 5 may reuse operand faces in a result.

## 2026-10-19: zero work excludes the filter

Pieces: 1, 2
Change: clarify that `max_work = 0` permits no decision work, including the f64 filter; a sign returns `NotExact` with a bound from its cached enclosure. Work accounting, defaults, public signatures and implementation behavior are unchanged.
Why: the architect ruled that piece 2's `not_exact_curve_and_point_distance_bounds` test incorrectly expected filter proofs with no work, despite piece 1 charging one work unit for an operation on doubles.
Migration: assertions expecting a filter-proven nonzero coordinate use `max_work = 1`, the smallest sufficient budget; assertions expecting zero-work uncertainty retain `max_work = 0`. Expected outcomes and decision counts stay unchanged.

## 2026-10-19: piece 4 deepened, with closed forms and branch-aware queries

Pieces: 4, 5, 6, 7
Change: `intersect_curves` and `CurveCurve` are removed; piece 5 finds where two curves on one face meet as the meeting of one curve with the other's second surface. Added: `CurveForm` and `curve_form`, `on_curve` (branch-aware membership), `same_curve`, `arc_point` (an exact point strictly inside an arc) and `TracePiece::taylor` (enclosures of Taylor coefficients). `tangent_at` now takes a `Decider` and returns `Result<[Vector3; 2], Refusal>`: `[arriving, leaving]` tangents, which differ only at the singular end of an intersection curve, where recognising the singular point takes decisions. An intersection curve built by hand on a pair of equal cylinders with meeting axes is refused with `intersect.unsupported-curve`. A curve-surface hit kept on a branch after a `NotExact` is stored as a meeting of the three surfaces, so no fact puts it on the curve. The refusal codes are now `intersect.unsupported-pair` and `intersect.unsupported-curve` only. The description fixes every configuration of plane, cylinder and sphere pairs, the representation rule (an analytic curve when its fields are rational, an intersection curve with a stated witness otherwise), the natural direction and order of every component, the construction facts recorded, and that `intersect_surfaces` decides with exact rational arithmetic, makes no `Decider` call and is never marked.
Why: blind tests need one exact answer per call; lines and circles with irrational offsets are common in coaxial and parallel poses and must still be exported and counted as lines and circles; nodes need one-sided tangents; measurement needs high-order enclosures for validated quadrature.
Migration: the piece 5 creator uses `intersect_curve_surface` and `on_curve` instead of `intersect_curves`, and reads both entries of `tangent_at` at vertices, passing its `Decider`. The piece 6 and 7 creators use `curve_form` and `trace`.

## 2026-10-19: piece 5 deepened, with section edges and rational point queries

Pieces: 5, cli
Change: `IntersectionGraph::curves` is replaced by `edges() -> &[SectionEdge]` (curve, end points and the faces of each body the arc lies on), and `points()` now returns `&[PointId]` of isolated meeting points only. `classify_point` takes `&mut Model` and `&RPoint3`: a rational query point lets every ray and probe be rational, so piece 4 can intersect them exactly, and the model may gain auxiliary geometry. New `face_point` gives a rational point strictly inside a face. `boolean.foreign-body` is removed (a handle from another model is a programming error, as `ARCHITECTURE.md` says); `boolean.imprint-ambiguous`, `boolean.merge-ambiguous` and `boolean.no-rational-point` are added. A `NotExact` whose degenerate branch would need a fact no proof gives (two points as one, a point on a curve, two curves or surfaces as one) is refused instead of marked, as DESIGN.md R5 requires. Passing one `BodyId` as both operands is a programming error. The description fixes the stage contracts, the classification rule imprint must satisfy, the touching and maximality rules and the cases that define correct behaviour.
Why: the stage outputs must be testable, and classification of curved faces needs exact probes, which rational data allows.
Migration: the front end passes no point to `classify_point`. Later operations reuse the stages through these accessors.

## 2026-10-19: piece 6 deepened

Pieces: 6, cli
Change: no signature change. The description fixes the accuracy rule (width at most `max(relative * m, absolute)` with `m` the smaller end magnitude), the exactness rule for polyhedra, this piece's own deterministic work count under the decider's policy, cones and tori as measured surfaces, and the new refusal codes `measure.centroid-not-separated` and `measure.budget-exhausted`. Measurement answers carry only `measure.*` marks; the body's own marks stay on the body.
Why: blind tests need a checkable width rule and stated outcomes when the budget runs out.
Migration: the front end merges the body's marks and the measurement's marks for the manifest.

## 2026-10-19: piece 7 deepened, with a fixed file layout

Pieces: 7, cli
Change: no signature change. The description fixes the header, the twenty fixed entities, the entity mapping, the depth-first numbering, the real-number text, the entities added for ring edges and loopless spheres and tori, the writing of voids and branches, and the three parts of `max_deviation`. `write_step` writes nothing on a refusal.
Why: byte-identical output and blind tests need one layout; readers need every face bounded.
Migration: none for the front end beyond reporting `added_vertices` and `added_edges` from the report.

## 2026-10-19: judge contract clarifications for edge kinds, self-touching points and slit edges

Pieces: judge interface, 5, cli
Change: `edges_by_kind` counts an edge by the shape of its curve, so a meeting of two surfaces that is a line, circle or ellipse counts as that kind, and `intersection` counts only branches with no closed form. One connected piece of material that touches itself at a point where its edges meet (the node of two cylinders tangent from inside, or of Viviani's curve) is refused with `boolean.self-touching-shell`, stage `sew`, like a piece touching itself along a curve. ND6 gains a maximality rule: no edge has both its coedges in one face, so a tangent segment or line left inside a face is removed with any vertex only it uses. A ring edge keeps its vertex when the vertex is a singular point of its curve. A merge of faces or edges that cannot be decided is refused with `boolean.merge-ambiguous` instead of merged with a mark. The example manifest's `max_deviation` is `0`, as piece 7's rule gives for an answer written with exact decimals. The file formats are unchanged, so the format versions stay at 2.
Why: lines and circles with irrational offsets are stored as intersection curves but are lines and circles; at such a node the two sides of the contact must share a face, so neither the shared-vertex nor the split-vertex form passes the topology check. Without the slit rule, the counts for a cylinder tangent to a face along a segment, or along a line from ring to ring, did not follow from the rules, and both readings passed `check`.
Migration: the front end counts edges by piece 4's `curve_form`. The judge may see this refusal in a few tangent poses of cylinder pairs and cylinder-sphere pairs.

## 2026-10-19: which operations see the deadline

Pieces: all
Change: `ARCHITECTURE.md` now states that only the `Decider` calls that do work (`sign`, `compare`, `enclose`, `real_roots`, `certify_root`) see the deadline, so an operation that reaches its answer without one answers after the deadline. Tests after the deadline accept the `numbers.deadline` refusal or the default-budget answer, unless a description fixes which.
Why: the `Decider` offers no way to read the deadline, so the earlier wording of pieces 4 to 7 promised refusals that a correct build could not give.
Migration: none in code. Test writers follow the rule in their deadline cases.

2026-10-19: Pieces 4 and 7, case corrections to piece 4 case 27's inner-tangency near misses and piece 7 case 10's void-face same_sense; both are case corrections with no rule or interface change.

## 2026-10-20: abandon work in progress at the answer deadline

Pieces: 1, 4, 5, 6, 7, cli, judge interface
Change: decisions check the deadline during budgeted work, not only on entry. `Decider::check_deadline()` provides a checkpoint without changing decision statistics. `TracePiece::{evaluate, derivative, taylor}` now take `&Decider` and return `Result<_, Refusal>` so series evaluation and subdivision can stop. `Ball::{pi_checked, atan_checked}` share the unbudgeted algorithms but check the deadline during their series. Measurement and STEP fitting propagate these refusals without retaining partial results. The judge limits paragraph now promises a `numbers.deadline` manifest and exit 0 within 1 s of the deadline, including work in progress.
Why: release cylinder-sphere questions exceeded a 2 s deadline by 4.751 s, and some exceeded a 50 s deadline without a manifest at 75 s. Checking only the start of a decision did not protect the 60 s outer limit. This supersedes the 2026-10-19 entry about which operations see the deadline.
Migration: pass the current decider to trace enclosures and propagate failures. Long loops outside numeric decisions can call `check_deadline()`. Use the checked pi and arctangent enclosures inside deadline-bearing operations. File formats, deterministic budgets and the determinism exception for deadline refusals do not change.

## 2026-10-20: an unaffordable exact attempt does not end a decision

Pieces: 1 (and every piece that decides through it)
Change: piece 1's rule 4 now says that an exact attempt (annihilating polynomial, separation bound, algebraic identity) whose proof would pass `max_work` is skipped while refinement continues. A decision ends `NotExact` only when the next refinement step would pass a limit, so a non-zero value whose sign refinement within both limits can show is always `Proven`. No signature, code or budget default changes.
Why: on main the sign loop stopped with `numbers.budget-exhausted` at 64 to 128 bits when a separation-bound attempt was too expensive. Piece 5's case 57 query at scale `10^-30` (case 13 bite) was refused although its value, about `10^-104` against coordinates near `10^5`, needs about 400 bits. With the attempt skipped instead, the query answers `Proven(Inside)` in 1.5 s.
Migration: piece 1 keeps refining after a failed or unaffordable exact attempt. Callers need no change; some `NotExact` outcomes become `Proven`.

## 2026-10-20: root isolation narrows quadratically

Pieces: 1 (and every piece that calls `real_roots`)
Change: `real_roots` narrows isolating intervals, including the narrowing that decides whether a root is rational, by quadratically convergent steps sized relative to the root, kept only when exact signs confirm the bracket, with bisection only as a fallback. Piece 1's case 40 gains a quartic with a 2,800-bit leading coefficient that must be isolated under the default budget. No signature, outcome, code or budget default changes.
Why: on main, when the rational candidates cannot be enumerated, each root is bisected down to width `1 / (4 q^2)`, with `q` the leading integer coefficient. For piece 5 case 23 at scale `10^-30`, a probe-plane quartic with `q` of 1,428 bits needed about 2,760 bisections per root, exhausted the 10^9 work budget and refused the boolean. A scale-relative safeguarded Newton narrowing isolated the same quartic, rationality check included, with 1.3e7 work, and the boolean answered with the expected counts.
Migration: piece 1 replaces bisection-only narrowing in `algebraic::isolate`. Callers need no change.

## 2026-10-20: real roots of degree up to 4 carry irreducible defining polynomials

Pieces: 1 (and every piece that calls `real_roots`)
Change: a root `real_roots` returns that is not rational carries a defining polynomial irreducible over Q whenever its square-free factor has degree at most 4, so a quartic that splits into two rational quadratics is split. Piece 1's case 40 gains such a quartic and an exactly zero expression on its roots that must be `ProvenCoincident` under the default budget. No signature, outcome, code or budget default changes.
Why: piece 4's hits of the case 29 bite's branch with the plane `z = 0` are roots of `(x^2 - a)(x^2 - b)` with rational `a`, `b`. Kept as roots of the quartic, an exactly zero square-root argument at those hits could not be proven within 10^9 work units, so `intersect_curve_surface` refused `numbers.sqrt-not-separated`, and so did piece 5's `classify_point` of a case 55 point with `z = 0`. With the quartic split into its quadratic factors, the same intersection answered with its two crossing hits using 7.8e6 work units.
Migration: piece 1 factors square-free factors of degree 4 with no rational root into rational quadratics where they exist (through the resolvent cubic, for example); rational roots are already removed. Callers need no change.

## 2026-10-20: quadratic towers, enclosure at a precision, and traces enclosed on demand

Pieces: 1, 4 (and pieces 5, 6, 7 through them)
Change: piece 1's rule 3 adds exact arithmetic in a tower of real quadratic extensions, tried before any annihilating polynomial or separation bound, for expressions whose leaves are rationals, square roots and roots of degree 2. Piece 1 adds `Decider::enclose_at(&self, x, precision) -> Result<Ball, Refusal>`, which encloses a defined `Real` at a working precision without changing decision statistics. Piece 4's trace pieces enclose their non-rational constants on demand at each call's precision through it, instead of at `max_precision_bits` when the trace is made. Piece 1's case 39 gains three tower cases.
Why: 198 of 232 refusals in the edge hunt were `numbers.deadline` at imprint. In its five smallest reproducers, single sign decisions on exact zeros in `Q(sqrt(5), sqrt(671))` or nested square roots took 13 to 24 s each, and region signatures re-froze trace constants to 16,384 bits at every projection and precision. With a quadratic tower, on-demand enclosure (approximated by a 512-bit cap) and one trace per edge, all five booleans answered unmarked in 5 to 31 s.
Migration: piece 1 adds the tower path and `enclose_at`. Piece 4 freezes trace constants lazily. Piece 5 traces each boundary edge once per region-signature call. Callers of `trace` need no change.

## 2026-10-20: exact periodic chart lifts and bounded geometric searches

Pieces: 1, 4, 5
Change: piece 1 embeds constructible quartic roots into the same exact quadratic tower, through a biquadratic depressed polynomial or a positive rational root of its cubic resolvent. Candidates must satisfy the original polynomial and its exact isolating interval. Piece 5 propagates exact integer periods when aligning cylindrical trace boxes. Direction and chart searches use 64 attempts, and inward face probes use 256 dyadic depths, independently of arithmetic precision. Piece 4 reuses denominator reciprocals in Taylor recurrences and skips arithmetic on exact-zero interval coefficients. No existing public signature, refusal code or `BudgetPolicy` default changes.
Why: independent quartic leaves prevented exact imprint zero proofs, while deriving a period shift from rounded midpoints introduced non-periodic errors and could make adjacent boxes disjoint at every precision. Precision-derived geometric loop counts magnified unresolved cases. Full-question profiles then placed most measurement time in high-order trace jets, including repeated reciprocals and exact-zero coefficient arithmetic.
Migration: callers need no change. Every successful numeric and geometric decision still has an exact proof; exhausting a geometric search refuses instead of choosing an unresolved result.

## 2026-10-20: shared validated trace jets

Pieces: 4, 6
Change: trace clones share a bounded recent cache of Taylor enclosures keyed by exact parameter midpoint, radius, working precision and order. Constant enclosures are cached separately for each requested precision. Measurement traces each edge once per scalar-property call, and shares the traces across volume, area and centroid within `mass_properties`. Each property retains its fresh arithmetic budget.
Why: full-question profiles placed most measurement time in repeated trace Taylor evaluation. Reusing an already certified enclosure removes that computation without changing the enclosure or its proof. Precision-specific constant caching also prevents call order from changing an enclosure.
Migration: no public signature, requested accuracy, refusal code or budget default changes. Cached requests still check the deadline. Cache eviction changes cost only.

## 2026-10-20: zero radicands retain working precision

Piece: 1
Change: `enclose_at` tries a bounded exact tower proof when a square-root argument's ball overlaps zero. A proved zero is enclosed as zero. Other known-defined arguments still use the square root of the ball's nonnegative part, without annihilator or separation-bound work.
Why: clipping an uncertain enclosure of an exact zero loses half the working bits at each nested square root. In the smoke probe, a nested zero root at 128 bits had width `3.1798e-10`; at 256 bits it still had width `5.6398e-20`. The exact proof gives width zero at both precisions and changes no decision statistics.
Migration: no signature, accuracy, refusal code or budget default changes. Positive radicands must never be clamped to zero solely because their balls overlap zero.

## 2026-10-21: shared canonical arc partitions

Piece: 5
Change: imprint splits section arcs and operand boundaries with the same canonical cut pool before identifying shared edges. A graph section can therefore map to several common topology edges; their union is exactly the original section.
Why: splitting an operand edge without splitting a coincident section left overlapping arcs with different endpoints. Their duplicate outgoing germs had no strict radial order. Applying one partition removes the duplicate arcs instead of adding a tie-breaker. This resolves the two overlap-driven imprint-order refusals; the remaining fourth-order contact requires the derivative extension below.
Migration: no public Rust signature changes. Consumers of `Imprinted::common_edges` must use the complete common-edge set, not assume a one-to-one correspondence with graph sections. Exact tangent and curvature predicates, genuine ambiguity refusals and default budgets are unchanged.

## 2026-10-21: binary ball storage and exact-zero queries

Pieces: 1, 4, 6
Change: piece 1 adds `Ball::is_exact_zero(&self) -> bool`, true only for the singleton `[0, 0]`, without allocating field copies. Ball midpoints and radii use private odd integer mantissas with binary exponents. Reciprocal rounding uses integer quotients of dyadics instead of normalizing rational endpoint fractions. Rational accessors, rounding, enclosure formulas and work charges are unchanged.
Why: validated Taylor arithmetic spent most of its time allocating and normalizing rational binary denominators, including copies made only to test a coefficient for zero. Binary exponents and exact quotient rounding remove those denominators and their repeated GCDs.
Migration: trace and measurement coefficient checks use `is_exact_zero`; the allocating helper functions are removed. Public rational accessors retain their signatures and canonical values. A before/after arithmetic probe retains identical midpoint/radius hashes, and the GMP/MPFR enclosure suites check containment independently.

## 2026-10-21: separation-conditioned root narrowing

Piece: 1
Change: root isolation uses the exact quadratic discriminant when rational divisor enumeration is incomplete. Irrational closed-form brackets are reserved for branches closer than the magnitude-scaled Newton guard; well-separated branches retain the ordinary isolator's brackets. Rejected magnitude-based Newton proposals can use interval Newton. Cross-factor separation skips only refinement targets that leave every interval unchanged, retaining the original power-of-two target grid. Exact signs still certify every bracket; the guard selects an exact algorithm, not geometric coincidence.
Why: translated clustered roots defeated magnitude-only proposals. The quadratic smoke changes from default-budget refusal at 999,914,807 work units to two proven irrational roots with 900,624; the analogous cubic changes from refusal at 999,941,208 to one proven irrational root with 214,834,275.
Migration: no public signatures, default budgets, rationality-width targets or irreducibility guarantees change. Independent GMP certificates check the root lists and brackets, including the translated rational cubic.

## 2026-10-21: regular higher-order radial contacts

Piece: 5
Change: after tangent and curvature ties, regular germs may use exact implicit Taylor coefficients through order eight. Both use the same transverse linear parameter, so unequal tangent speeds cannot change the order. Singular or unresolved germs and shared components still refuse.
Why: the remaining contact case has a common sphere and internally tangent cylinders of radii three and two. With transverse parameter `z`, both branches have `y = z^2/4 + z^4/64 + ...`, but their `x-3` coefficients are `-z^4/96` and `-z^4/64`. Second-order predicates cannot distinguish them.
Migration: no public Rust signature or default budget changes. The internal radial sorter receives the model and each germ's source. Tests derive third- and fourth-order conic coefficients independently, vary tangent speed, reverse comparison order and retain refusal for duplicate germs.

## 2026-10-21: retained exact probe incidence

Piece: 5
Change: probe-plane constructions record the incidence they establish, and rational face probes retain the exact side-query proof. For an unmarked surface pair with one curve, existing incidence on both surfaces proves query incidence on that curve. Pairs with multiple components retain the branch-membership predicate. Construction provenance separates auxiliary probe geometry from intrinsic seed supports, preventing later probes from recursively following earlier probe sections while retaining the facts themselves.
Why: rebuilding an already-known conic in a different algebraic coordinate representation forced another zero elimination and exhausted the decision. The missing information was the existing incidence proof, not a geometric tolerance.
Migration: no public API or budget changes. Derived facts retain exact provenance. The two-generator meridian regression compares public curve membership with independent line-distance predicates, so support membership cannot be mistaken for membership of the wrong component.

## 2026-10-21: bounded exact root-certificate reuse

Piece: 1
Change: each decider lazily holds at most 128 nonlinear root certificates in a FIFO memo. A hit requires exact equality of every rational polynomial coefficient. Complete roots, multiplicities and isolating intervals are reused unchanged; hash collisions never identify polynomials. Shared polynomial keys avoid copying their coefficients into both lookup and eviction indexes. Reuse charges the original certificate work and retains deadline checks.
Why: a full native cylinder case profile spent 35.75 seconds in 3,286 root calls for 63 polynomial fingerprints; 34.78 seconds repeated successful fingerprints. Re-isolating the same exact polynomial dominated the remaining deadline rather than needing more precision or a larger budget.
Migration: no public signatures, default budgets, repeated-decision work guarantees or geometric tolerances change. The memo is bounded and decider-local. Independent GMP certificates cover a working set larger than the memo, and a warm-cache test retains work identity and expired-deadline refusal.

## 2026-10-21: enclose certified root intervals without re-proving their endpoints

Piece: 1
Change: algebraic enclosure checks the certified interval width before evaluating the defining polynomial at an endpoint. That sign is computed only when bisection is required. Refining calls keep the same bisection sequence; already-narrow calls return the same midpoint/radius bounds without an unused polynomial evaluation.
Why: the residual native profiles spent 48/99 and 25/99 samples in algebraic enclosure. Deep isolating intervals can have thousands of denominator bits even when the caller asks for a modest working precision. Re-evaluating their large defining polynomials consumed time and could spend an enclosure budget without improving its bounds.
Migration: no signatures, default budgets or rounding rules change. The regression first isolates `sqrt(2)/(2^607-1)`, then reaches width `1e-220` under a separate 10,000-work policy; independent GMP signed-square inequalities certify the returned bounds.

## 2026-10-21: shared algebraic enclosure certificates

Piece: 1
Change: immutable irrational-root snapshots share their polynomial, interval and a bounded memo of eight working-precision enclosures. A cache hit returns the exact ball produced by the original bisection at that precision and replays its work cost. An insufficient allowance takes ordinary computation to preserve partial-work accounting. Refinement and defining-polynomial replacement do not inherit another snapshot's memo.
Why: adjacency coefficient checks repeatedly enclosed the same algebraic leaves. Public root-list memoization did not retain those bisection certificates. Snapshot sharing also avoids copying large polynomial coefficients and interval endpoints whenever an algebraic value is cloned.
Migration: no public signature, default limit, rounding rule or algebraic-equality contract changes. No higher-precision substitution occurs. Independent GMP checks cover budget-limited cold/warm results and a cached request still refuses an expired deadline.

## 2026-10-21: certify stalled algebraic-zero decimals

Piece: 7
Change: after decimal refinement stalls on an enclosure containing zero, an exact piece 1 coincidence proof may establish the decimal `0.` with zero error. A containing interval alone never permits this result.
Why: a nonnegative enclosure of an algebraic zero can have a positive midpoint at every precision. Relative decimal refinement then stalls even though the defining polynomial proves the value exactly zero.
Migration: no signature, rounding rule or default budget changes. Previously successful certificates take the same path. Regressions use a cubic algebraic identity and a tiny positive perturbation to distinguish exact zero from a small unresolved value.


## 2026-10-19: what a boolean stage record's statistics mean

Pieces: 5, cli
Change: in `StageRecord::stats`, every counter other than `max_bits` is what the stage added to the `Decider`'s statistics, and `max_bits` is the `Decider`'s running `max_bits` when the stage ended, an upper bound on the stage's own peak. Piece 5 reads statistics only through `stats()`, never `take_stats()`.
Why: piece 1's `Decider` keeps only a running maximum precision and offers no way to restore cleared statistics, so a stage's own `max_bits` cannot be observed without clearing the caller's totals.
Migration: piece 5 snapshots `dec.stats()` around each stage. Its test writer checks case 53 as now written. No signature changes.

## 2026-10-19: STEP export refuses when the budget allows no bound

Pieces: 7
Change: new refusal code `step.budget-exhausted`. The export refuses with it, writing nothing, when the budget allows no proven bound on a value or spline fit the file needs, or no enclosure narrow enough to write a number within one unit in its last written place. A spline fit with a proven bound that misses `curve_tolerance` is still written with the mark `step.tolerance-not-reached`.
Why: the description gave no outcome for an export whose budget reaches no bound at all, and the export must not invent a `numbers.*` refusal. Piece 6 already has the matching `measure.budget-exhausted`.
Migration: the piece 7 creator returns the code in those cases. The piece 7 test writer accepts it in case 18 and its zero-work variant.

## 2026-10-19: a step may name one body twice

Pieces: judge interface, cli
Change: a step's `args` may name the same operand or step twice. `X ∪ X` and `X ∩ X` are `X`, and `X − X` is empty. `catallax-cli` answers such a step by that set rule, without calling piece 5's `boolean`, whenever both arguments resolve to the same `BodyId`. The result carries `X`'s marks once, and the step makes no decision. The question format is unchanged, so its version stays at 2.
Why: the question format already accepted repeated names, but piece 5 treats one `BodyId` passed as both operands as a programming error, so such a valid question made the program panic.
Migration: the CLI builder applies the set rule. Pieces 1 to 7 are unchanged.

## 2026-10-22: publish Catallax and retain the version 2 exchange

Pieces: all, cli, exchange
Change: workspace packages, Rust crate paths and directories use `catallax-` and `catallax_`. `catallax-cli` builds `catallax` and `exact-brep` from the same main. Benchmark and dev-only deadline environment identifiers use `CATALLAX_`. Protocol strings `exact-brep-question/2`, `exact-brep-manifest/2` and `kernel.name = exact-brep` remain unchanged.
Why: the B-rep implementation becomes the starting code for a new public MIT kernel with one recipe and two evaluators. Rename packages without breaking the existing exchange.
Migration: Rust callers use the new crate names. Existing exchange clients may keep invoking `exact-brep`. A future protocol rename requires a versioned contract change. The new refusal rule, 1 nm necessity bound and operation-derived persistent names are documented requirements with explicit implementation gaps, not geometry changes in this publication.

## 2026-10-22: recipe format `catallax-recipe/1` and crate `catallax-recipe`

Pieces: recipe (new); the judge interface is unchanged
Change: `docs/RECIPE.md` defines the recipe language, version 1: operation definitions, origin-based face and edge names, references, rebuild failures and the canonical JSON file `catallax-recipe/1`. The new crate `catallax-recipe` holds its data model: `Recipe`, `Step`, `Operation` and their parameter types, `Origin`, `Role`, `EdgeName`, `VertexName`, `OriginPattern`, `FaceRef`, `EdgeRef`, `EdgeSelection`, `Number`, `RecipeError` with codes `recipe.<reason>`, `Recipe::parse`, `Recipe::check` and `Recipe::to_text`. It depends only on piece 1, and nothing depends on it yet.
Why: PRD section 10 step 3. Both kernels and SolidDraft need the same definitions and names before either kernel implements an operation.
Migration: none. No piece, the CLI or the exchange contract changes. The question format stays the judge's input.

## 2026-10-22: analytic fields and stock solids in recipe version 2

Pieces: recipe, field; B-rep pieces and the exchange interface are unchanged
Change: add `catallax-field` with `Plan`, a separate `Face`/`Origin` registry, `Workspace`, `Batch`, tentative `Approximate` results, and proven `Query`, `Location`, `Value`, `DistanceBound`, `Direction` and `FaceQuery` results. Recipe `Operation` gains exact `Block`, `Cylinder`, `Sphere`, `Cone` and `Torus` variants. Their face names use the existing start/end and side roles with fixed segment ids. `VertexName` accepts one origin for singular vertices such as cone apices, but still rejects empty or repeated origins. The schema identifier becomes `catallax-recipe/2`.
Why: both kernels need the same stock-solid definitions before comparisons can cover what the B-rep already builds. Residuals, true distance enclosures, conservative geometric bounds and unresolved Boolean zeros must remain distinguishable.
Migration: update operation matches for the five variants. Explicitly relabel reviewed version 1 recipes as version 2, then parse/check/write them; existing operation meanings and exact numbers are unchanged. Version 1 is not accepted as an alias. Persist recipe origins rather than `FaceId` or tape indices, and establish a registry face's surviving domain with `face_query`. Handle non-exhaustive query enums, especially `Indeterminate`; Boolean distance upper bounds are unbounded unless the point is certified on the boundary. No existing B-rep caller, CLI question, manifest or STEP format changes.

## 2026-10-22: algebraic geometry definitions and conic integration order

Pieces: 1 to 7, recipe evaluation, field and display consumers; documentation decision, implementation pending
Change: approve existing `Real` for defining frames, positions, dimensions and poses; add `RealPoly` beside rational `MPoly` and certified rational-degree sine/cosine construction. Analytic constructors validate through the caller's decider; chart methods accept that context and return enclosure failures explicitly. Definitions containing `Real` lose structural equality/hash traits; geometry comparisons use proofs. Add `EllipticCylinder`, preserving actual sweep trims. First land a separate conic preparation change with `Real` dimensions, parabola/hyperbola `CurveForm` variants and `curve_is_closed`; then widen all directions and consumers in one shared foundation. Generic quadrics and derived sweep/revolution surfaces remain deferred with explicit coverage obligations.
Why: legal torus sections, open cone sections, 45-degree revolutions, irrational-radius arcs and oblique extrusions exceed rational definitions. Field predicates read the same surfaces, so changing only the B-rep constructors would split the geometry contract. A surface-only type change would also leave rational polynomial, equality, metric and export assumptions unsound.
Migration: the approved sections in the piece documents supersede their retained baseline signatures where stated. Preserve rational input conveniences and successful rational exchange bytes in the conic preparation. The foundation adapts existing intersection, topology, boolean, measurement, STEP and field consumers together. Feature lanes rebase once onto that merged base; each then adds certified coverage and independent regressions. Read `ARCHITECTURE.md` for lane order and alternatives. No code, recipe schema, CLI contract, default budget or 1 nm policy changes in this documentation PR. `docs/JUDGE-INTERFACE.md` remains unchanged.

## 2026-10-22: B-rep display meshes and OBJ export

Pieces: display mesh, consuming 1 to 4; the judge interface is unchanged
Change: new crate `catallax-mesh` exports `mesh_body`, `MeshOptions`, `Mesh`, `Triangle`, `FaceMesh` and `EdgeMesh`. A mesh has shared binary64 positions, face-tagged triangles, per-face proven chord bounds and per-edge certified polylines. `Mesh::write_obj` writes shared vertices and `face_<id>` groups. Options specify an exact positive millimetre chord tolerance and a deterministic vertex limit. Invalid options return `mesh.invalid-options`; resource, chart-resolution or coordinate-precision failures return `mesh.precision-or-budget`. Numeric trace/deadline errors propagate and input marks remain attached.
Why: PRD section 6.5 requires meshes from the B-rep for display and later background comparison with the field kernel. Display approximations must not replace the exact body or claim the 1 nm agreement guarantee.
Migration: new library consumers may mesh an existing `Model` and `BodyId`, then pick by `Triangle::face`. Arena face IDs are not recipe origin names. No existing caller, STEP byte format, question format or manifest format changes.

## 2026-10-22: exact conic dimensions and curve closure

Pieces: 4, 5, 7, cli
Change: `CurveForm::Circle.radius_sq` and `CurveForm::Ellipse.major_sq` and `minor_sq` are `Real` rather than `Rational`. `CurveForm` adds oriented `Parabola` and single-branch `Hyperbola` forms. Piece 4 exports `curve_is_closed(geom, curve, dec) -> Result<bool, Refusal>` and boolean callers use it to distinguish ordered open curves from closed loops.
Why: even a horizontal torus section can have irrational squared circle radii. Cone sections also need open conic forms, so curve kind cannot stand in for closure.
Migration: construct exact dimensions with `Real::from` and compare them through `Decider`, retaining exact rational arithmetic when both expressions expose rational values. Primitive geometry definitions remain rational. Analytic storage requires rational fields; otherwise retain the intersection definition. STEP retains the existing rational square-root path and accepts exact real conic dimensions. Parabola and hyperbola forms specify the next pair implementations; this change adds no supported surface pair. The CLI schema and existing curve-kind names remain unchanged. Its `intersection` edge tag still means no recognised conic form; exposing open conics requires a separate versioned exchange decision.

## 2026-10-22: named B-rep sketch features and recipe version 3

Pieces: recipe, topology, boolean, features, measurement; the judge exchange is unchanged

Change: add `catallax-features` with exact sketch validation, odd-nested profiles, signed extrusion, analytic revolution, recipe stock-solid inputs and explicit booleans. `Operation::ExtrudeThrough` defines a two-sided through-all cut whose finite construction caps are proved outside the target. The schema becomes `catallax-recipe/3`. `Model::add_named_face`, `face_origins`, `edge_name` and `vertex_name` expose canonical recipe provenance without changing `Face` layout or unnamed legacy builders. Boolean splits retain origin sets; maximal face merges union them.

Builds retain an immutable recipe snapshot so operation definitions remain available after analytic reduction. Full revolutions reject isolated axis contacts that would pinch a shell; valid partial sectors remain buildable. Straight cone generators ending at an apex contribute exact zero to the boundary integrals, avoiding a removable singularity without raising measurement budgets.

`catallax_boolean::make_maximal` reuses the sewing merger on constructed boundaries while preserving their proven shell/lump nesting and singular loops. Half-turn caps joined along the axis merge with both origins; artificial axis edges and degree-two vertices disappear. The existing boolean stages retain their previous classification and refusal behaviour.

Why: recipes need actual B-rep features and persistent operation names. Rational analytic representations are an explicit intermediate implementation boundary, not new impossible-input rules. RECIPE.md section 7.14 lists the remaining legal-input bugs and their specific unsupported-geometry codes.

Migration: explicitly relabel reviewed version 1 or 2 recipes as version 3, then parse/check/write. Earlier operation meanings are unchanged. Match the new operation variant; field evaluation of sketch features and through-all remains unsupported. Use `catallax_features::evaluate` for named B-rep recipe builds, then the existing measurement and STEP APIs on its returned model and body. No command, question format, manifest format, legacy unnamed geometry or STEP format changes.

## 2026-10-22: rational closed-form cone and torus boolean paths

Pieces: 4, 5, 6, 7; the exchange contract is unchanged
Change: add plane–cone conics and generator rays, circular plane–torus sections, coaxial conic pairs and certified sphere/tube separation. `supports` now reports an implemented path for all five current surface kinds; it does not promise every pose. Curve membership and equality retain nappes and ray/branch domains, open conics have finite traces, and odd higher-order torus contacts remain crossings. Boolean charts and probes handle conical faces, apex loops and complete circular torus regions. Measurement integrates torus parallels by their exact Fourier mean and cone generators by their proven zero boundary form.
Why: countersinks, apex cuts and circular torus cuts need complete intersect/imprint/classify/keep/sew/measure/STEP paths, not blanket kind refusals or quadrature of removable singularities.
Migration: signatures, stored rational frames, CLI kind tags, budgets and STEP schema are unchanged. Missing general traces return `intersect.traced-pair-required`; irrational direction storage returns `intersect.algebraic-direction-required`; unimplemented noncircular torus charts return `boolean.torus-chart-required`. These are temporary implementation gaps. Open conics are not relabelled as `intersection` in the legacy manifest; their exchange exposure still requires the separately versioned contract decision. Rebase onto the coordinated algebraic foundation before widening these definitions.
Unresolved nappe and new conic component selectors refuse hit construction with the originating numeric code, rather than retaining the wrong component with a mark. The existing bounded branch marks for legacy plans remain unchanged.

## 2026-10-22: field profile sweeps and through-all cuts

Pieces: field, features, recipe
Change: `catallax_field::Plan::compile` accepts sketch profiles, signed `extrude`, two-sided `extrude_through` and rational-frame full/partial `revolve`. Sketch validation uses `catallax_features::Profile::new`; feature supports and trims are evaluated without calling a B-rep builder. Source face candidates retain segment aliases, split domains, reversed cut orientations and 180-degree merged-cap names.
Why: both recipe evaluators must build the same feature definitions and report their creating origins.
Migration: no Rust signature, recipe format or command-line protocol changes. Sweep values use `Value::Residual` and `Approximate::true_distance_formula = false`; unsigned distance upper bounds remain unbounded off the certified boundary. Callers must not treat a sweep residual as true Euclidean distance. The specific unfinished `features.unsupported-*` codes from RECIPE.md section 7.14 now apply to field compilation as well.

## 2026-10-22: sampled two-kernel display comparison

Pieces: 5, field queries, recipe consumers, display comparison
Change: add `catallax-compare` with `Controller`, revision tokens, scoped `Options`, complete `Report` values and explicit agreement/disagreement/unresolved outcomes. Add piece 5 query helpers `face_contains` and `surface_witness`, and `Plan::is_proven_empty` in the independent field evaluator. No existing geometry builder, recipe schema, question, manifest or STEP signature changes.
Why: the background display check needs certified distances in both directions, exact trimmed-face associations, face origin sets and cancellation that cannot publish stale or partial work.
Migration: existing callers are unchanged. Comparison callers begin a revision per edit and dispatch its check on their own worker. Read `Controller::latest` only as a revision-tagged snapshot. Treat regional and sampled agreement as scoped display evidence, never as a 1 nm or whole-boundary certificate. `surface_witness` supplies an exact rational chart projection, not a closest point; measure its distance and call `face_contains` before using it as an upper-distance witness.

## 2026-10-22: comparison consumes field profile sweeps

Pieces: field queries, display comparison
Change: add `Plan::visit_feature_candidates` as a query-only, cooperatively stoppable source of rational supporting-face candidates for compiled sketch sweeps. Comparison now accepts sketch, extrusion, through-all and revolution steps in the common rational support subset. Candidate generation propagates conservative field-result bounds through the expression tree to sample through-all walls without finite construction caps.
Why: the field evaluator gained these operations while the comparison library was being built. They must use the same two-way checks as stock solids rather than retaining obsolete unsupported-operation reports.
Migration: call `face_query` before accepting a candidate as a surviving boundary point. The visitor returns `Ok(false)` when its callback stops, not a complete-coverage claim. `field.sample-options`, `field.sample-unbounded`, `field.sample-role` and `field.sample-algebraic-radius` are sampling failures, not geometry refusals or agreement. Existing builders, query meanings, recipe schema and exchange formats are unchanged.

## 2026-10-22: composable geometric error decision

Pieces: 1 to 7, recipe evaluation, field, mesh, compare
Change: [ARCHITECTURE.md](ARCHITECTURE.md#composable-geometric-error-decision-of-2026-10-22) defines a whole-result geometric certificate, typed local conversions, operation transfer proofs and a strict bound below 1 nm. It reserves `numbers.geometric-gap-unproven` and `numbers.geometric-gap-exceeded` for implementation.
Why: inherited mark logs record local uncertainty without proving its effect over the final solid, a recipe chain or an exported representation.
Migration: implement the decision in the documented piece order on the shared algebraic foundation. Audit every producer and consumer before certifying marked results; keep measurement and display uncertainty separate. This is a documentation contract, with no runtime API or code change. The command-line contract, legacy manifest meanings and work-budget policy remain unchanged. A later implementation must log its concrete API migration; certificate exposure through the CLI needs a separate versioned exchange decision.

## 2026-10-22: recipe version 4 definition contract

Pieces: recipe specification and implementation plan. Documentation only; no crate, CLI, exchange format, work-budget or STEP interface changes in this commit.

Change: RECIPE.md specifies `catallax-recipe/4`. Loft now uses rational arc charts, uniform section stations and common parameter refinement for unequal segment counts. Thicken uses local projection sectors at sharp selected joins and regularised normal sweeps with straight boundary rims. Paths admit stationary round turns, explicitly directed reversals, helices and rational splines. Sweep exposes interior-generated boundary patches with `envelope.path` origins. Solved algebraic sketch coordinates/radii retain rational polynomial coefficients and exact root selectors. Through-all gains one-sided tools; `extrude_to` and `extrude_next` define bounded first-hit extrusions and their cap/curtain names. Shell regularisation and the empty-kept-face case are explicit. Fillet/chamfer roll over neighbouring faces by default, with explicit keep-neighbours setbacks and consumed-face history. Draft retains neutral-connected material with kept-region overrides. `scan_import` is a non-solid reference value with separate reconstruction deviation reports.

Why: both kernels need the same sets, branches, names and refusal rules before feature implementation. Rational loft correspondence changes an existing operation's meaning, so version 3 must not silently acquire it. Exact algebraic coefficients do not by themselves provide general derived carriers or certify transcendental paths.

Migration: the current crate and evaluators remain version 3 and must reject version 4 until its schema lands. A later implementation must add the new variants, fields, root construction validation and origin roles, with explicit version conversion. Review arc-bearing lofts and sharp-join thicken edits for changed geometry; retain version 3 when the old definition is intended. Unchanged operations may migrate after reference validation. Canonical root definitions normalise their polynomial, not the geometry to a rounded coordinate. Existing code and fixtures are unchanged here.

Feature patterns expand to saved ordinary steps; outward and two-sided shell/thicken constructions compose existing operations. Import replacement rewrites the saved hash, decisions and references in one edit. Rebuild never refreshes `near`, and authored splits retain the directed-start segment id. These editorial rules allocate no hidden kernel names. The selected product policies are roll-over with explicit neighbour preservation, automatic neutral-connected draft with region overrides, and reference-only open scans. General variable-radius laws and nonplanar parting-line operations remain separate future schema extensions. No implementation coverage is claimed by this specification bump.
## 2026-10-22: deterministic verification exchange

Pieces: numeric accounting and CLI verification
Change: add `catallax-verify stream`, `catallax-request/1`, `catallax-response/1` and `catallax-manifest/3`. Question, recipe/3 and two-kernel comparison sources share a deterministic whole-operation allowance. Add `with_operation_budget` for synchronous request scopes. Manifest/3 defines distinct elliptic-cylinder, parabola and hyperbola tags.
Why: the legacy 50-second deadline lets machine load choose between an answer and a refusal. Verification needs deterministic outcomes and separate latency evidence.
Migration: existing `catallax` and `exact-brep` commands, question/2 and manifest/2 remain unchanged. Adopt the new binary explicitly. Question adapters parse but ignore legacy `deadline_ms`; `max_work` is the new allowance. Do not dispatch numeric work across threads inside a request scope. Parallelize independent requests. See [VERIFY-INTERFACE.md](VERIFY-INTERFACE.md) for streaming, cache identity and cancellation behavior.

## 2026-10-22: shared algebraic foundation

Pieces: 1 to 7, field, display mesh, two-kernel comparison and cli consumers; judge formats and kind tags unchanged
Change: implement the algebraic-definition decision across the existing consumers. `RealPoly` stores ordered `Real` coefficients, with certified degree, equality, resultants and univariate root selection; rational polynomials retain the existing `MPoly` backend. `sin_cos_degrees` constructs certified algebraic sine/cosine values from exact rational degrees. `Decider::rational_value` optionally proves an exact rational value without replacing its expression definition; `None` is inconclusive, not a proof of irrationality or definedness.
Why: frames, poses, analytic dimensions and their support equations must describe the same exact geometry in both kernels, including irrational rigid motions and oblique circular extrusion.
Migration: use `Point3`, `Vector3` and `Real` for stored definitions, including `PointDef::Explicit` and `CurveForm` directions. Analytic constructors and `Pose::new` take a final `&mut Decider`; `Rotation::from_real_quaternion` adds checked algebraic rotation input. The rational quaternion constructor remains available. `matrix()` exposes `Real` entries; `rational_matrix()` and `rational_components()` borrow exact rational entries when available. Definitions containing `Real` have no mathematical equality/hash traits; use certified geometry comparisons. Construction identity can prove cache reuse only.

`Pose::apply_rational_point` and `apply_rational_vector` preserve native rational-input arithmetic without temporary expression leaves. They return `None` when the stored frame does not expose rational entries; they do not round or decide whether an algebraic expression happens to be rational. The existing recipe feature constructors use this path; runtime field placements and ray samples remain algebraic.

`Pose::from_rational` accepts a validated rotation and a native rational translation without a numeric proof. Rational recipe sampling uses this constructor and the native point transform.

Charts take the caller's decision context and return an outer `Result`, retaining their former inner optional result where applicable. Implicits and curve equations return `RealPoly`; consumers use `rational()` or `into_rational()` for the borrowed or consuming rational fast path, and certified operations otherwise. No rational conversion rounds an algebraic definition.

`EllipticCylinder` adds its normalized transverse implicit, periodic/linear chart, plane intersections, measurement over the actual oblique trims, and STEP `SURFACE_OF_LINEAR_EXTRUSION` over an `ELLIPSE`. Equal semi-axes reduce to a circular cylinder before construction. New-kind display meshes and new recipe coverage remain separate work.

Field plans retain algebraic support polynomials, gradients, ray germs and rigid placements. `compile_checked`, `proven_checked`, `face_query_checked` and `transformed` accept the shared decider. Convenience queries use the default policy and return an explicitly indeterminate, unbounded result when a numeric proof fails; checked queries preserve the numeric refusal. Rigid transforms preserve the existing metric formulas; no affine-stretch distance shortcut is introduced.

The comparison consumer uses `RealPoly` support equations and gradients, certifies the first nonzero normalisation coefficient, and compares supports through the shared decider. Rational equations retain the original `MPoly` arithmetic. Rational-only witness and feature-candidate APIs do not round algebraic placements: an unavailable native candidate frame returns `field.sample-algebraic-placement`, a sampling failure rather than a geometry refusal or an agreement claim.

Algebraic field coefficient, metric and ray-germ enclosures use `Decider::enclose_budgeted_at`, preserving the caller's precision and work limits. The existing `enclose_at` remains the explicitly unbudgeted, known-defined working helper. Native rational enclosures keep their original arithmetic path.

Intersection-plan caches retain input definition nodes. Their keys therefore use `Surface::shares_definition_nodes`; value-equal independent definitions cannot inherit another query's expression-DAG aliases and change charged proof work. `Real::same_node` is a storage-identity query, not a mathematical equality test.

This supersedes the rational-only primitive-storage statement in the preceding conic-dimensions entry. Question/manifest version 2, supported CLI kind tags, numeric limits and deadlines are unchanged. Elliptic cylinders and open conics are not exposed through new or remapped CLI tags.

## 2026-10-22: algebraic definitions in the verification exchange

Pieces: CLI verification and algebraic geometry consumers
Change: the manifest/3 publisher recognises `elliptic-cylinder` under the preceding deterministic-verification decision. Parabolas and hyperbolas retain their distinct manifest/3 tags. The shared-algebraic-foundation entry's CLI restriction applies to the legacy question/2 and manifest/2 exchange, whose tags and refusals remain unchanged.
Why: the verification exchange landed while its geometry consumers were being widened; its explicit versioned vocabulary must not inherit the legacy publication restriction.
Migration: select extended publication only for `catallax-verify`. Keep legacy commands on their existing vocabulary. Request-scoped work accounting also covers the new checked algebraic proofs; neither exchange's numeric limits or default allowance changes.

## 2026-10-22: retain curve-meeting and root-refinement proofs

Pieces: 1, 5
Change: explicit curve meetings name both proven carriers in their point definitions. Line substitution retains an actual common-polynomial certificate: rational equations use their GCD, and a sole nonzero algebraic equation certifies its roots directly. Otherwise shared predicates check every carrier equation. Analytic conic candidates retain their exact equation checks without repeating curve membership. A finite meeting with either support exhausts the candidates for an intersection branch; the other support and component predicate remain mandatory. Certified branch membership retains its original proof as incidence. Root enclosure certificates retain exact bisection endpoints and the lower-end sign, so higher-precision calls can continue the same sequence when their sign-filter precision agrees, replaying the cold-evaluation work.
Why: later boundary queries discarded already-proven curve incidence, and increasing enclosure precision repeated exact polynomial evaluation over the same bisection prefix.
Migration: no public signatures, default budgets, deadline or exchange-format changes. A completed refinement retains the same bisection sequence and rounding. Cold and warmed caches retain identical work charges under the selected sign algorithm. Refined snapshots and changed defining polynomials retain separate certificates.

## 2026-10-22: skip zero and repeated Taylor products

Piece: 4
Change: interval Taylor multiplication and division skip coefficient ranges proved identically zero. Squared expressions and square-root recurrences compute symmetric products once, while retaining the original order of rounded additions and subtractions.
Why: high-order traces used by face membership and measurement repeatedly multiplied the same coefficient pairs or scanned zero tails.
Migration: no public signatures, trace definitions, orders, precision or rounding rules change. A before/after public-API probe retained identical midpoint, radius and precision fingerprints for 57,456 coefficients across circles and quartic branches, Taylor orders zero through twenty and four working precisions.

## 2026-10-22: filter algebraic bisection signs

Piece: 1
Change: when the odd part of a bisection denominator exceeds 512 bits, root enclosure tries a budgeted polynomial ball at the caller's working precision, capped at 512 bits, before exact rational Horner evaluation. Smaller odd denominators retain the direct exact path. Only a strict interval sign selects a side; zero-containing enclosures retain exact evaluation. Bisection points, root selection and final rounding stay unchanged.
Why: a certified bracket can have very large odd denominators without needing equally large working precision. Exact polynomial values at every midpoint exhausted the budget even for the real root of `x^3 - 2`.
Migration: no public signature, default allowance or deadline changes. A regression brackets `cbrt(2)` with 4,096-bit odd denominators and requires width `1e-30` within one million work units; GMP cubic inequalities check containment. The cold/warm precision-history matrix checks work identity and insufficient-allowance fallback.

## 2026-10-22: schedule separated point coordinates first

Piece: 2
Change: point equality first schedules coordinate differences with a cheap nonzero rational value or a cached conservative interval excluding zero. If these do not separate, metered enclosures at the existing interval widths precede exact-zero comparisons. Enclosures only schedule queries; the decider still certifies each conclusion under its original policy and deadline. Coincident points retain their original coordinate order, and the distance bound remains the symmetric Euclidean bound.
Why: canonical-point searches attempted expensive algebraic coincidence in an earlier coordinate even when another coordinate already separated the points.
Migration: no public signatures, geometric definitions or default policies change. Predicate regressions retain their incidence, bounded uncertainty, domain and deadline checks; fixed query-count assertions that prescribed coordinate order are removed.

## 2026-10-22: seed polynomial content from a small coefficient

Piece: 1
Change: primitive integer normalization starts its coefficient GCD from a nonzero coefficient of minimum bit length, stops when the content is one, and avoids dividing an already primitive positive-leading polynomial by one.
Why: annihilator construction spent time computing content from very large constant terms even when a unit coefficient already certified primitive content.
Migration: no public signature, policy or normalized coefficient changes. Existing polynomial and root suites retain their rational-content, sign, multiplicity and independent-oracle checks. A 65,536-bit monic-polynomial smoke probe retains every coefficient.

## 2026-10-22: certify polynomial coprimality in a finite field

Piece: 1
Change: polynomial GCD first tries a fixed 31-bit prime. It concludes coprimality over the rationals only when every coefficient denominator is invertible, both degrees are preserved, and finite-field Euclid ends at a nonzero constant. Otherwise the existing exact polynomial GCD runs unchanged.
Why: square-free annihilators repeatedly entered large integer remainder and content calculations merely to establish that the GCD with their derivative was one.
Migration: no public signatures, coefficient domains, default budgets or deadline changes. Normalized GCDs stay unchanged. Regressions cover a common factor lost by modular degree reduction, noninvertible rational denominators and coprime rational polynomials whose modular images share a factor.

## 2026-10-22: retain singleton-zero enclosure certificates

Piece: 1
Change: a known-defined immutable zero singleton, or a zero singleton retained after checked refinement, proves coincidence without algebraic elimination or a nonzero-root separation bound. The existing algebraic proof category is used.
Why: sign decisions discarded exact zero enclosures and could exhaust their allowance reconstructing a zero already certified by interval arithmetic.
Migration: no public signatures, default allowances or deadline changes. A ball merely containing zero proves nothing. Definedness must be established before using a singleton, including when zero multiplies an expression with a hidden undefined divisor. The regression fails before this change at one work unit; empty allowances still stop, and tiny nonzero enclosures remain non-coincident.

## 2026-10-22: retain contact proofs across shared supports

Piece: 5
Change: adjacency probes retain zero constant terms proved by shared source equations or typed point/curve-to-surface incidence facts, even when another boundary curve uses a different identity. A whole carrier equation and its exact negation both certify the same zero set.
Why: expanding a known contact into real polynomial coefficients discarded its construction proof.
Migration: no curve-membership conclusion is inferred for another selected component. Whole-path exclusion and component checks remain mandatory; no public signature or geometric definition changes.

## 2026-10-22: try separation certificates before contact proofs

Piece: 5
Change: adjacency probes try boundary equations with zero-excluding cached coefficients first. Endpoint checks likewise try coordinates whose source and endpoint enclosures are disjoint before potentially coincident coordinates.
Why: an exact contact in the first coordinate could trigger expensive algebraic elimination even when another coordinate immediately separated the whole path.
Migration: cached bounds select only the order of existential checks. The same exact whole-path and selected-component predicates still decide adjacency; no public signature, default allowance or deadline changes.

## 2026-10-22: share zero-system query scheduling without reordering bounds

Piece: 2
Change: curve membership now uses the same separation-first scheduling as point equality. A fixed-size pending array retains each equation's original index, so reordered decisions still feed the original equation units into distance bounds.
Why: a coincident first curve equation could enter expensive elimination before another equation could establish separation.
Migration: no public signatures or default policy changes. Cached signs and metered enclosures select the query order only; the decider still certifies every conclusion. Existing analytic-curve, intersection-carrier, zero-work bound, undefined-coordinate and deadline tests pass.

## 2026-10-22: prune and prioritize polynomial critical pairs

Piece: 4
Change: Buchberger reduction uses the chain criterion only for pairs with previously established standard representations. Remaining pairs are selected by leading-monomial LCM degree, then monomial and generator indices.
Why: redundant FIFO reductions caused coefficient swell in posed finite curve/surface intersections.
Migration: coefficient operations and domains are unchanged; the criterion uses only leading monomials and resolved-pair identities. Reduced monic bases, root ordering and stored geometry definitions retain the same contracts. A cyclic three-equation regression checks both exact real solutions under all six generator orders; the complete scale and placement family passes.

## 2026-10-22: scale preliminary coefficient enclosures

Piece: 5
Change: the probe-coefficient sign filter first tries the cached conservative interval width when it is looser than the refinement target. Further refinement is bounded relative to that width as well as by the existing absolute target.
Why: at large geometric scales, an absolute width alone forced unnecessary precision before testing a coefficient's sign.
Migration: this changes only a preliminary metered enclosure request. A zero-containing ball still proves nothing, and the exact fallback is unchanged. No geometric tolerance, default allowance, public signature or deadline changes.

## 2026-10-22: omit Taylor coefficients unused by panel certificates

Piece: 6
Change: integration programs lazily compile center and remainder dependency plans into fixed-size boxed slices, keeping both plans independently lazy without spare-capacity metadata. Symmetric center integration needs even output coefficients; its whole-range certificate needs only the order-20 output coefficient. Products, variable quotients and roots retain every lower coefficient needed by their recurrences. Constant divisors retain their zeroth coefficient and requested numerator coefficients.
Compiled node, output-index and exact-constant definitions are immutable boxed slices too. The builder counts retained nodes while marking dependencies and allocates their final buffer once; compaction adds no allocation or coefficient arithmetic. Prepared planar area state retains only its numerator and squared normal. Optional open-torus program bundles are indirect, so their payload is not embedded in every curved face.
Why: evaluating every output coefficient forced arithmetic that never entered the integral or its remainder bound.
Migration: requested coefficient arithmetic keeps its order and outward rounding; definedness checks and nonsymmetric angular integrals remain intact. No public signature, requested accuracy, default allowance or deadline changes. Regressions use exact antiderivatives and keep undefined zeroth-order operations visible.

## 2026-10-22: consume temporary dyadic magnitudes

Piece: 1
Change: ball products and midpoint-rounding residuals consume fresh dyadic values when taking absolute values. Midpoint/radius normalization also retains owned, already-precision-sized results instead of copying their integer mantissas. Precision reduction of borrowed balls keeps a borrowed path.
Why: radius formulas and interval normalization copied temporary results before immediately discarding them. An odd normalized mantissa changes whenever significant bits are discarded, so the representable-center case needs no rounded copy or error subtraction.
Migration: public signatures, arithmetic, rounding and work charges are unchanged. A 57,456-coefficient public jet workload retains identical midpoints, radii and precisions; independent GMP/MPFR tests check containment.

## 2026-10-22: share immutable evaluation-memo enclosures

Piece: 1
Change: one real-expression evaluation shares immutable ball values between memo entries and their callers. The memo is dropped before returning the owned root ball. Square-root clipping replaces only its local value.
Why: memo insertions and hits copied arbitrary-precision midpoint and radius mantissas.
Migration: memo scope, expression traversal, domain checks, rounding and work charges are unchanged. No public ownership or thread-safety contract changes. The complete numbers suite passes, including cold/warm budget boundaries and known-defined square-root clipping.

## 2026-10-22: avoid rational differences in diamond-chart comparisons

Piece: 5
Change: cylindrical projection compares each diamond phase against shifted anchor thresholds instead of normalizing a rational difference. Corner bounds are computed once, and four phases use a fixed array.
Why: chart subdivision repeatedly normalized fractions only to compare their order.
Migration: phase quotients, corner order, strict wrap cutoffs and outward hull rounding are unchanged. A before/after smoke retains identical results for 2,916 scaled boxes at 8, 64, 128 and 256 bits. Permanent checks independently bound all four quadrants, a periodic seam and an axis, and reject a box containing the projection pole.

## 2026-10-22: restrict straight contact quadrics using their source proof

Piece: 5
Change: a verified straight planar probe with proven source contact uses `f(t) = t*c1 + t²*(f(q)-c1)`, where `c1 = grad(f)(q)·p + 2*f(q) - grad(f)(q)·q`. The coefficients are affine in the algebraic source point.
Why: expanding quadratic source-coordinate products reconstructed a zero already supplied by the contact proof.
Migration: only exact carrier checks and existing incidence proofs enable the reduction. Contact polynomials and endpoints use `RealPoly` and `Point3`; the shared decider certifies degree and endpoint contact without rational-only coefficient assumptions. Other paths retain their original expansion, and whole-path separation and component checks remain mandatory. Regressions independently derive simple, double and contained contacts, reject extraneous algebraic carrier candidates, and retain non-rational quadric contact/interior roots. A standalone smoke checks nine coefficients and fifteen exact path values.

## 2026-10-22: retain owned multivariate accumulators

Piece: 1
Change: owned polynomial addition retains its left term tree and moves the right terms into it. Entry-based accumulation moves vacant coefficients directly, computes occupied sums once, and removes exact cancellations.
Why: curve substitution repeatedly cloned growing term trees and copied incoming arbitrary-precision coefficients before discarding them.
Migration: operator signatures, coefficient domains, term order and normalized results are unchanged. Borrowed addition retains its independent result. Targeted independent polynomial tests cover substitution, cancellation, derivatives, resultants and certified systems; no precision or work-policy changes.

## 2026-10-22: retain connected-exterior intersection proofs

Piece: 5
Change: an exact, nonsingular intersection of bounded operands without void shells retains the proof that every result shell is outer. Each operand has a connected exterior; those exteriors meet outside the operands' bounds, so their union is connected too. Sewing no longer reconstructs that fact through rational region sampling and ray-root isolation.
Why: tiny, rigidly placed cylinder/sphere intersections exhausted root-isolation work while rediscovering a shell role already implied by the operation and input topology.
Migration: selected face senses and topology construction are unchanged. Marked classifications, singular contacts and void-bearing operands retain geometric shell certification. No public signature, coefficient domain, default allowance or deadline changes. The complete 456-Boolean default scale/placement family passes with all original topology and material witnesses, and runs normally again.

## 2026-10-22: coplanar booleans and same-body set rules

Pieces: boolean engine, field and two-kernel comparison
Change: the boolean engine accepts the same body as both operands. Subtraction returns empty; union and intersection normalise the retained material through sew. The field classifies complete local plane arrangements, including cancelled common supports, edges and corners. Opposite coincident primitive halfspaces certify empty intersections during compilation. Comparison no longer excludes repeated operands.
Why: flush unions and coplanar subtraction left removed field supports unresolved, and recipe subtraction from the same body reached a distinct-body assertion.
Migration: no format or public signature changes. Empty identities still have no boundary samples; separating halfspaces now supply another exact emptiness proof. General curved degeneracies and unbounded Boolean upper distances can remain unresolved. The legacy question exchange is unchanged.

## 2026-10-22: question construction without measurement

Pieces: CLI library, application callers, measurement
Change: `catallax_cli::build_question(&Question, &mut Decider)` returns a `BuiltQuestion` with its owning `Model`, result `body`, construction `marks` and boolean `per_step` statistics. It performs no property integration, topology summary or export. `BuildError` retains the refusal and optional boolean step and stage.
Why: an application needs a usable body before requesting proven volume, area or centroid. A later measurement refusal must not discard that body.
Migration: library consumers that need a body use `build_question`, then call `catallax_measure::volume`, `area` or `mass_properties` on demand. Recipe applications keep using the already unmeasured `catallax_features::evaluate` path. The legacy commands, answer functions, verification defaults, requested enclosures and STEP contract are unchanged. Verification questions can continue opting out with an empty `measure.properties` list.

## 2026-10-22: SolidDraft conversion and application build entry

Pieces: recipe, features, mesh, application
Change: new `catallax-app` crate exposes `soliddraft::convert`, `build` and `build_soliddraft`. Results carry the exact `Built`, marked display mesh with face and edge picking handles, complete face-origin sets and optional persistent edge/vertex names. SolidDraft conversion returns the original-tag to canonical-step table. Typed errors separate invalid source, not-yet conversion or operations, step-local build refusals, topology faults and display refusals; `Error::code()` returns the stable code. No legacy CLI or exchange format changes.
Why: SolidDraft needs one entry from recipes to named, displayable bodies without inventing topology identities or taking mesh geometry as the model. Fillet, chamfer, shell and unresolved conversion constructions remain explicit not-yet results.
Migration: SolidDraft's backend consumes this entry in its own repository. Map returned handles explicitly to backend ids, preserve merged origin sets, and persist names rather than handle indices. `docs/RECIPE.md` section 6.3 describes the output contract and current conversion limits. Recipe booleans also apply exact repeated-operand set identities before invoking the distinct-operand intersection graph, so a body combined with itself no longer panics.

## 2026-10-22: reuse failed sphere-gauge enclosure attempts

Piece: 6
Change: grouped volume, area and centroid measurement retain request-local sphere-gauge domain failures by face, subdivision depth, rectangle/Taylor mode, precision and gauge. A later property skips an attempt only when all six gauges already failed the identical area-gauge denominator check and every boundary arc uses the same trace integration path. Property-specific analytic arcs and second moments do not use this cache.
Why: centroid measurement repeated the same failed domain enclosures that volume and area had already resolved by subdivision. Those failed attempts consumed its measurement allowance before the first finite centroid enclosure.
Migration: no public signature, refusal code, requested accuracy, default allowance or deadline changes. Successful panel arithmetic and retry ordering are unchanged. Individual failed gauges are not skipped inside a successful attempt, so its quadrature-dominance flag remains unchanged. A failed enclosure is never treated as proof that a geometric point or region is absent.

## 2026-10-22: proven centroidal inertia on request

Pieces: measurement, application and question measurement publication
Change: `MassProperties` includes `inertia: Option<InertiaTensor>`, a symmetric 3-by-3 enclosure matrix about the centroid at unit density. Empty bodies have no centroid or tensor. `inertia` requests only the tensor; `centroid` requests only the centre. `volume_area_centroid` preserves the original grouped measurement without second-moment work. `ApplicationBody::mass_properties` and `ApplicationBody::inertia` borrow the built result, so a refusal does not discard it.
Why: applications need proven inertia for mass-property reporting and equivalent B-rep bodies, rather than mesh estimates or guessed values.
Migration: callers of `mass_properties` receive the additional measurement and at most one additional width mark. Call `volume_area_centroid` to retain the original property set and budgets. Questions opt in with `measure.properties: ["inertia"]`; the matrix uses model axes and length units to the fifth power, with negative products of inertia off the diagonal. Default question measurements, other verification answer formats and STEP definitions are unchanged.

## 2026-10-22: whole-face definition deviation certificates

Pieces: field, verification CLI, verification driver
Change: `catallax_field::measurement` exposes `CentrePath`, `FaceDefinition`, `MeasurementOptions`, `FaceMeasurement` and `measure_face`. It proves lower and upper bounds on the worst deviation of a B-rep face from its own rolling-path or offset-support definition. Exact faces need no record. The recipe verification adapter accepts an opt-in `catallax-face-check/1` wrapper and reports `catallax-face-measurements/1`, including per-face local tolerances, first-mesh timing and resolved origins after a declared scalar edit. Field `face_measure` log entries carry those local certificates instead of null.
Why: fillet prototypes and later tolerance consumers need whole-face, two-sided definition checks. A Boolean residual, boundary point query or display chord bound cannot substitute for them.
Migration: existing requests, recipe/3, legacy manifests and STEP exports are unchanged. Opt-in callers attach the exact definition to named result faces and consume the attained interval and `converged` flag. Supporting-chart covers can conservatively include trimmed-away regions. Path arcs retain their exact branch and endpoint definitions. Invalid controls and unsupported covers use `field.face-measurement`; metadata and reference errors use `cli.face-definition`, with numeric and geometric refusals propagated. Determinism comparisons also exclude the diagnostic `face_measurements.display.first_mesh_ms`, while retaining certificates, origins and work counts. No fillet or offset construction, tolerance cap or merge policy is added.

## 2026-10-23: planar display trims retain sub-ulp gaps

Pieces: display mesh
Change: planar faces now triangulate their trims with certified `Real` orientation and containment predicates before rounding display positions. Hole bridges and ears retain every sampled boundary vertex and shared edge index. Curved faces keep the periodic CDT path.
Why: a valid flange union leaves a gap smaller than one binary64 ulp between a cap's outer loop and a hole. Projecting rounded positions first made those constraints overlap and split without orientation labels. The exact body and Boolean are unchanged.
Migration: no API, model, STEP or chord-tolerance changes. A planar triangle's bound is the maximum certified endpoint displacement, since its exact support is affine. Sub-ulp strips may have coincident rounded positions but retain distinct indices and closed incidence. The minimized two-block union and multiple-bore regressions cover these boundaries.

## 2026-10-23: exact application queries and optional display bounds

Pieces: measurement, Boolean point queries, meshing and application
Change: `catallax-measure` adds `face_centroid`, `edge_length`, `edge_midpoint`,
`bounding_box`, `face_bounding_box` and `BoundingBox`. `ApplicationBody` borrows its
existing build for each query. Bounds select analytic extrema with exact trim
membership and use budgeted second-derivative enclosures for general branches.
`face_contains` accepts exact `Point3` through a source-compatible generic argument;
its rational caller and classifier are unchanged. Measurement now depends on the
Boolean crate without introducing a dependency cycle.

Display options add `angular_deviation`, `max_edge_length` and an optional
`max_aspect_ratio` request. Facet normal certificates include analytic normal
variation across each parameter triangle and coordinate rounding. Length checks
use the exact rational values of returned binary64 positions. Face reports give
an aspect-ratio bound and the requested certified angular bound. Unmet aspect
requests appear in `Mesh::unavailable_controls`, with the face and reason.
Exact planar trims retain their incidence through shared dyadic edge refinement.
`TracePiece::evaluate_exact` retains the algebraic positions of boundary samples
before display rounding, including samples introduced by length controls.
The aspect bound is `None` when rounding collapses a facet to zero area; requested
aspect and angular controls then appear as unavailable rather than finite bounds.
Planar triangulation proposes indexed CDT connectivity, then certifies the exact
convex hull, every triangle orientation and the oriented trim-boundary incidence.
Collapsed or uncertified proposals retain the exact trim clipper. Outward-rounded
binary64 interval filters accept only proven signs; inconclusive predicates use
`Real` arithmetic.
Sine/cosine range reduction uses a stored 4096-bit Machin enclosure of pi with
outward truncation and a series fallback above that precision. The constant has
no numeric work charge or per-request cache state. Direct pi queries are unchanged.
Edge subdivision reuses its certified endpoint enclosures down the recursion.

Why: applications need exact geometric summaries and explicit display guarantees,
not quantities inferred from triangles or controls silently discarded.
Migration: mesh option literals should use `..Default::default()` for optional
controls. An edge midpoint is halfway along arc length, including closed edges
with their start vertex or carrier seam. Bounds contain no display chord padding;
an empty body has no box. Existing mass-property requests, question defaults,
answer formats and STEP definitions are unchanged.
