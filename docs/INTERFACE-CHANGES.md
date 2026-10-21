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
