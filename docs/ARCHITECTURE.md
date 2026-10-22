# Catallax architecture

Catallax's master model is one recipe evaluated by an exact B-rep kernel and a signed distance field kernel. Each operation shares one mathematical definition. [PRD.md](PRD.md) defines that direction; [DESIGN.md](DESIGN.md) governs the B-rep implementation.

The seven pieces below implement the B-rep half, exported without source history into this public repository. [RECIPE.md](RECIPE.md) defines the shared recipe language. `catallax-recipe` holds its data model, and `catallax-field` evaluates stock solids and booleans with separate named face records. `catallax-compare` checks their sampled display boundaries in both directions; conversions and the finer pre-export check remain unimplemented. [FIELD.md](FIELD.md) defines the field query boundary; [PLAN.md](PLAN.md) orders the remaining work.

The target is reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels. The former 1-in-10,000 programme and grid-cell product gates are retired. Mathematical regression tests remain; speed per click matters for the interactive CAD tool.

## Exact geometry definitions: decision of 2026-10-22

Adopt algebraic analytic definitions. Use the existing `Real` expression type for defining positions, directions, dimensions and rigid motions. It already contains only algebraic operations and certified algebraic leaves. Keep rational input types and rational arithmetic as the fast path. This is an approved implementation contract, not a claim that the code below already implements it. The migration sections in pieces 1 to 7 take precedence over their retained rational baseline descriptions.

| Decision | Reason | Alternative and cost |
| --- | --- | --- |
| Use `Real`, with proven constructor invariants | Retains construction expressions and shared roots without eagerly eliminating each coordinate | Store every field as `Algebraic`. Equality still needs proof, and eager conversion can grow polynomial degrees and discard useful construction structure |
| Add `RealPoly` beside rational `MPoly` | Algebraic frames produce algebraic implicit coefficients; intersections and field signs must read the same definitions | Keep rational polynomials by eliminating coefficients immediately. This introduces conjugate solutions and requires branch constraints at every use |
| Add `EllipticCylinder` now; retain the five existing analytic surface kinds | Oblique circle extrusion needs it. Current revolve profiles use lines and circular arcs with an axis in the sketch plane | Add every quadric now. More kinds require chart, domain, singularity, intersection, integration and export contracts without a current recipe operation requiring them |
| Defer generic `Extruded` and `Revolved` surface variants until profiles need them | Lines and circles have analytic reductions that reuse exact predicates and distance methods | Add derived surfaces now. This retains a direct sweep parameter but duplicates analytic carriers and needs cross-kind equality and a general certified intersection path immediately |
| Retain the operation definition and names through analytic reduction | A support surface alone loses trims, source segments and branch choices | Store only the reduced carrier. This is smaller but cannot reconstruct the feature or distinguish surviving patches |

Other quadrics join when an operation requires them, with a chart/domain and proof contract, not as an unchecked polynomial escape hatch. Generic derived surfaces are required for spline profiles, sweeps and lofts whose supports do not reduce analytically. They must retain immutable profile/path references, exact laws, domains and orientation; their charts, predicates, intersection certificates, integration and export must land together. A point evaluator alone is insufficient. A valid patch of a non-ring circular revolution also needs a domain-aware carrier before that case can be claimed complete; the stock ring-torus constructor must not be used to declare that patch impossible.

The costs are larger expression graphs, algebraic equality decisions and potentially high-degree elimination. Rational degrees have algebraic sine and cosine, but their degree is unbounded as angle denominators grow. Retain rational and quadratic fast paths, shared immutable expressions and construction facts. Do not rationalise, snap a frame, raise global budgets, or eagerly compute a minimal polynomial for every field. Measure construction, predicates, intersections and repeated field queries separately. Eager canonical numbers would simplify cache keys at the cost of paying elimination even for display-only queries.

Determinism comes from fixed construction and monomial orders, certified root selectors, deterministic work accounting and ordered result assembly. Expression identity is not mathematical equality. Caches must not change proof outcomes or work allowances; thread scheduling must not choose roots, branches or topology. The separate verification exchange below resolves the wall-clock deadline conflict for callers that adopt it. Widening types does not resolve it, change the 1 nm rule or make unsupported geometry impossible.

### Integration order

Use one shared foundation, with the following merge order. Feature branches consume its final merged base once; they do not each introduce a competing geometry representation.

1. Land this documentation decision.
2. Land the separate conic preparation change first: `CurveForm` dimensions become `Real`, parabola and hyperbola forms and the closure query are added, and all current consumers compile. Existing successful rational exchange fixtures must retain byte-identical STEP and geometric manifest content, excluding elapsed-time fields. This step adds representation, not a claim of complete cone/torus intersections.
3. Land a single coordinated algebraic foundation across numbers, shapes and every existing consumer, including field. It adds `RealPoly`, widens frames/poses and stored analytic dimensions, adds elliptic cylinders, and supplies certified angle construction. It includes the consumer adaptations needed to preserve the existing rational subset. Publish this merged commit as the common base. Keeping the rational subset usable avoids a period in which either evaluator needs its own fork of the definitions.
4. Rebase the cone/torus, sketch/extrude/revolve, display-mesh, fix batch 5 and field feature lanes once onto that base. Independent fixes and the already-authorised rational recipe subset may land before step 3 if they do not widen interfaces; fold those commits into the common base. Later corrections must preserve the agreed interfaces or receive a new recorded decision.

| Lane | Adaptation after the common base |
| --- | --- |
| Cone/torus booleans | Use algebraic coefficients and dimensions; preserve nappe and component selectors; implement open-conic membership/order/finite traces; use proven closure instead of treating every non-line as a ring |
| Sketch, extrude, revolve | Build exact sweep vectors and rational-degree rotations; reduce supported sides to analytic definitions; retain trims and origin sets; remove rational-only refusals as each full operation path gains coverage |
| Display meshes | Consume certified charts, derivatives and finite trimmed traces for all supported kinds; cache rounded vertices only; retain face origins and deviation bounds; do not rediscover topology from samples |
| Fix batch 5 | Keep fixes on the shared predicates and closure query; remove rational-only coefficient assumptions in contact, probe and classification paths; retain independent regressions and construction facts |
| Field kernel | Migrate support polynomials, gradients, ray germs and pose comparisons in step 3; extend recipe coverage on the same definitions afterwards; certify metric bounds separately from polynomial signs |

Cone/torus and recipe coverage may then land independently, with consumer and agreement regressions for each added case. Merge their completed operation paths before declaring display/field agreement for those operations. Uncoordinated per-lane widening is the rejected alternative: it gives earlier isolated patches but forces repeated rebases and leaves field and B-rep interpreting different types. No lane may report a valid input as impossible because its implementation is unfinished. Keep specific temporary failure reasons and coverage tests until repaired. The command-line contract in `JUDGE-INTERFACE.md` remains unchanged. Its fixed kind tags cannot describe elliptic-cylinder faces or open conics. The separate manifest/3 exchange in [VERIFY-INTERFACE.md](VERIFY-INTERFACE.md) assigns those results distinct tags. No silent tag remapping is authorised.

## Composable geometric error: decision of 2026-10-22

Adopt one proof of geometric error for the whole result, relative to its exact recipe definition. A mark records why a decision stopped; it is not that proof. This is the implementation contract for [PLAN.md](PLAN.md)'s error-policy item. It supersedes the inherited rule below that automatically takes every `NotExact` as coincidence. This documentation decision changes no code, work-budget policy or [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md) contract.

| Decision | Reason | Alternative and cost |
| --- | --- | --- |
| Certify the represented solid against the exact recipe result | Bounds remain meaningful after operations, trims and changes of representation | Certify only individual predicates. Small residuals can select a distant branch or remove material |
| Require a whole-result bound strictly below `1/1_000_000` mm | Many individually small decisions can exhaust the same allowance | Give every operation a fresh 1 nm. A long recipe can drift without limit |
| Retain typed local evidence and a shared dependency graph | A consumer can verify scope, units, conditioning and inherited error | Sum or take the maximum of the mark log. Neither operation supplies missing geometric proofs |
| Require a transfer proof for each operation | Boolean boundaries and topology can change discontinuously near contact | Assume booleans are non-expansive because field expressions use `min` and `max`. Those expressions need not be true distance |
| Keep property and display approximation bounds separate from model error | Area, volume and display chords have different meanings | Treat every bound as millimetres. This accepts dimensionally invalid evidence |

### What the whole-result certificate proves

Let `R` be the regularised solid defined by the exact recipe prefix, including its exact branch and domain rules. Let `B` be the valid solid represented by the returned body. Neither refers to a mesh or a numerical sample. For a nonempty bounded regular closed solid `S`, define `s_S(x)` as the true Euclidean distance to its boundary, negative inside and positive outside, in millimetres. Define

```text
D(B, R) = sup over x in Euclidean 3-space of |s_B(x) - s_R(x)|
tau     = 1/1_000_000 mm
certificate: D(B, R) <= E, where E is a nonnegative exact rational and E < tau
```

This is an oriented geometric bound. It bounds two-way boundary distance and guarantees the same material classification wherever `|s_R(x)| > E`. Boundary Hausdorff distance alone would miss an inside/outside reversal. The supremum is a proof obligation, not an instruction to sample all space. An implementation may discharge it by analytic bounds or a complete certified cover of trimmed patches, cells and their exterior, with material orientation checked. A finite sample set never discharges it. Intermediate points, curves and open patches use their own scoped correspondence bounds; they are not solids to which `D` applies.

Two empty results have distance zero. An empty result and a nonempty result have no finite certificate under this policy. A bounded decision may not erase the last component and call the empty output close. Prove exact emptiness or refuse. Multiple bodies exported as an assembly need per-body certificates and a matching proof; overlapping bodies are not implicitly replaced by their union.

`E = 0` with a proof means exact geometry. A positive conservative enclosure remains bounded even if the actual result happens to be exact. An empty log, a zero cached residual or an exact-looking file is not an exactness proof. Query and output enclosures may still have width when model geometry is exact.

The certificate must identify the immutable recipe revision, represented body or query, units, covered domain, exact target, bound and proof dependencies. Preserve the original decision code, residual enclosure, affected entities and branch choice beside each geometric conversion. Use stable recipe/definition identities and ordered proof records, never pointer addresses or worker completion order. Shared dependencies form an acyclic graph. Geometry remains stored by definition; approximations and their certificates are separate, explicit representations. A mark never authorises a false construction fact or an invalid shell whose edge merely lies near its face.

### From a local mark to millimetres

Every geometry-changing use of `NotExact` needs both a conversion and a proof of the effect of taking the degenerate branch on the entire affected domain. Each accepted coincidence substitute must itself have a proved geometric gap strictly below `tau`, as well as preserving the whole-result cap. Later cancellation cannot excuse an over-limit local choice. A conversion at one witness does not cover a face, a trimmed curve or a body. Piece 1 cannot infer geometric units from a `Real` expression. The geometric caller owns the conversion and proves every denominator, domain restriction and conditioning bound used in it.

| Evidence | Required conversion and scope |
| --- | --- |
| Distance `b` | Use `delta = b` only for the objects and domain named by its certificate. Distance to a carrier plane is not distance to its trimmed face. Distance to two surfaces separately is not distance to their intersection. Establish a branch-preserving projection or another geometric correspondence before using it for either purpose |
| Angle `a` in radians | Prove the angle convention and a finite lever arm `R`. Rotation about a fixed pivot/axis moves every covered point by at most `R*a`; `2*R*sin(a/2)` may tighten this with a certified enclosure and a valid angle range. Translation uncertainty adds separately. Unbounded supports need finite trims. A normal-angle bound does not bound intersection displacement without a separate intersection-conditioning proof |
| Scalar parameter `b` | Declare its dimension and exact map `X(t)`. Prove `sup norm(dX/dt) <= L` on the whole parameter interval, then `delta <= L*b`. For several parameters, `sum L_i*b_i` is a conservative bound. This is pointwise correspondence evidence; topology and the whole-result claim still need proof |
| Raw quantity `b`, including polynomial values, dot products and determinants | Supply an explicit residual-to-geometry theorem. For `n dot (x-o)`, `b/m` bounds plane distance when `norm(n) >= m > 0`. For a proven squared point distance, an outward bound on `sqrt(b)` is a distance. Rescaling an equation must not change the certified geometric meaning. A generic polynomial residual cannot become a distance just by relabelling it |
| Measurement-only scalar or raw quantity | Keep the property's units and its enclosure. An area width is in mm², a volume width in mm³, a centroid coordinate width in mm. Such an evaluation does not move geometry and contributes zero additional model displacement. Its correctness and requested accuracy are separate obligations. If a property later controls geometry, that consumer must supply a sensitivity proof |

A local inverse bound, for example a residual divided by a derivative lower bound, requires existence of the intended root and validity throughout a certified neighbourhood. A nonzero gradient at one point does not suffice. Intersections need branch existence and an inverse/Jacobian bound over the chosen domain; near-parallel surfaces can amplify arbitrarily small residuals. Singular contacts need a specialised bound, such as a certified square-root modulus, or an exact decision. An angle or scalar with no finite domain or conversion remains unproven. No generic conversion may set its contribution to zero.

Observed enum terminology is `Distance`, `Angle` and `Value`, not four separate variants: [marks.rs](../crates/catallax-numbers/src/marks.rs), `BoundKind`. The policy distinguishes a scalar with known units from an untyped raw `Value`; that distinction requires new evidence, not a new spelling in the existing CLI manifest.

### Composition and the five boolean stages

An operation must prove a transfer relation for the exact operation `F`, the represented inputs `B_i`, their ideal inputs `R_i` and its output `B_out`:

```text
D(B_out, F(R_1, ..., R_n)) <= T_F(E_1, ..., E_n; input certificates) + delta_F
E_out = an outward rational upper bound of that right-hand side
```

`delta_F` covers new bounded choices relative to `F(B_1, ..., B_n)`. `T_F` covers inherited error through the operation. A direct certificate against `F(R_1, ..., R_n)` may replace this decomposition. No default `T_F = max(E_i)` or `sum(E_i)` is valid for arbitrary booleans. Even an exact operation on already bounded inputs needs a transfer proof. Near contact, a small input displacement can create or destroy a component or change which distant patch survives.

For certified consecutive representations of the same target, use the triangle inequality and add their errors. For a complete collection of corresponding patches, a maximum is permitted only after proving that the cover, material orientation and gluing establish the whole-solid metric. A Lipschitz transfer `sum L_i*E_i` is permitted only with proved constants on the whole input uncertainty domain. Otherwise refine or refuse. Certifying a smaller bound directly against the original target is allowed; discarding history is not.

| Stage | Required proof before its result can be used |
| --- | --- |
| Intersect | Retain every relevant meeting, component selector and finite trim. Convert local residuals to bounds for the intended meeting locus. Prove the effect of an omitted or merged meeting, including its possible downstream classification changes |
| Imprint | Prove that the split arrangement covers each input face, that cyclic ordering and incidence are valid, and that merged cuts stay within their certified correspondence. Carry intersection evidence forward. Counting it a second time does not create a new error |
| Classify | Certify material labels outside the allowed uncertainty region, with oriented local evidence inside it. Taking a tangent or on-boundary branch must not select an unbounded or distant face. An unresolved sample cannot label a whole face without a coverage proof |
| Keep | Derive the selected material from the operation's regularised set semantics for both operands. Account for lost and retained patches, including thin walls, cavities and possible emptiness. This stage must not treat a small local gap as proof that a component is dispensable |
| Sew, including maximal merging | Prove a valid represented solid, complete face/edge correspondence and orientation, then close the whole-result certificate against the exact recipe prefix. Account for every new merge or representation change. Publish the body only after `E_out < tau` is proved |

Intermediate stage objects must carry scoped certificates and unresolved proof obligations explicitly. A reusable stage API must not expose such an object as a certified body. Check any available composed bound at each stage and stop before returning an over-limit answer. Check again when closing the result. If a stage cannot bound its downstream effect yet, it may retain alternatives for later proof, but may not irrevocably discard one and publish it as certified.

Recipe evaluation follows the same rule at every body-producing node. Exact rigid motions preserve `E`; a certified uniform scale multiplies it by the absolute scale. Other transformations and features require their own transfer proof. Reusing a body in two operands shares evidence but does not remove its influence through either operand position. For an additive transfer proof, sum both contributions unless a joint proof establishes their dependence. Deduplicate proof identity, never numeric bounds or reason strings. Exact identities such as `A subtract A = empty` may discharge inherited error for that output while retaining the audit record. Separate equal-looking bodies require an equality proof before applying such identities.

There is no per-operation reset and no fixed allocation of 1 nm among stages. A chain with successive certified displacements of `0.4 nm` and `0.4 nm` can have `E = 0.8 nm`; another `0.2 nm` reaches the limit and fails the strict test. Cancellation in the displacement may reduce that bound only if proved. Sharing a subexpression in a recipe must not allow a later operation to overlook its earlier error.

### Refusals and publication

Reserve these stable codes in piece 1 so all consumers report the same failed invariant. The implementation entry must identify the producing piece, recipe step and boolean stage where applicable.

| Code | Trigger |
| --- | --- |
| `numbers.geometric-gap-unproven` | A needed local conversion, transfer theorem, valid representation, complete scope or whole-result certificate is missing or cannot be proved within the assigned work allowance |
| `numbers.geometric-gap-exceeded` | The best certified local substitution gap or whole-result upper bound available for acceptance is greater than or equal to `tau` |

The second code means the available proof cannot meet the cap; it does not claim that the true error exceeds 1 nm. Refine or take an exact branch within the assigned deterministic allowance if possible. If neither succeeds, return the specific refusal and no successful body, property hand-off or export. Never clamp the bound, increase the geometric cap, silently drop marks or label a bounded coincidence as `ProvenCoincident`. These are unfinished-kernel refusals under PRD section 3.2, not proof of impossible input. Existing domain and coverage refusals remain more specific when they occur before a bounded candidate exists.

The certificate is independent of the diagnostic log. Every piece must reject a legacy marked input without sufficient evidence, preserve certificates through passthrough paths, and validate imported or deserialised evidence before trusting it. A producer cannot self-certify merely by constructing a public mark. Final publication checks the whole certificate, including any representation change after the last boolean.

### Measurements, STEP and display

Piece 6 must state whether an enclosure describes `B` or the ideal recipe solid `R`. Shared mass properties must enclose `R`. Numerical quadrature and rounding error for `B` add to certified geometric uncertainty, not to the model's displacement ledger. From the signed-distance certificate, the symmetric difference of `B` and `R` lies within the `E`-tube of `B`'s boundary. A certified upper volume `V_tube` for that tube therefore bounds `|volume(B)-volume(R)|`. Use a complete tube bound, not the unsupported formula `area*E` near corners or contacts. Bounds on coordinate magnitude in the tube also bound first-moment differences. Divide moment intervals by a volume interval proved strictly positive to enclose the centroid. If zero cannot be excluded, preserve `measure.centroid-not-separated`.

A small distance bound alone gives no useful general surface-area bound: arbitrarily corrugated surfaces can be close and have very different areas. Area needs certified patch correspondence and Jacobian bounds, or direct integration over the exact recipe definition. If unavailable, refuse that property with `numbers.geometric-gap-unproven`; do not widen it by an invented multiple of `E`. A property enclosure that misses requested accuracy remains an unfinished accuracy failure even when its geometry certificate passes. Reporting such a diagnostic interval must not count as satisfying the hand-off.

Piece 7 must certify the solid described by the emitted numbers, fitted curves, analytic surfaces and trims. A local curve-fit bound alone is insufficient. For certified export displacement `E_export` relative to `B`, require

```text
E_file = E_model + E_export < tau
```

Include decimal rounding, fitted branches, endpoints, surface parameters, edge/face compatibility and any added seams in `E_export`. Maxima may combine independent covered entities only with a whole-solid correspondence proof; consecutive fitting and rounding add unless a joint bound is proved. An exact model does not give the file a free additional 1 nm. Requested export accuracy can be stricter; it cannot relax this cap. Build and check the certificate before publishing the file. Refusal must not leave a success report or a partial artifact presented as usable. Future export accuracy levels and seam support need their own exchange decisions.

Display meshing retains a separate mesh-to-`B` bound `E_mesh`, including position rounding and complete triangle coverage. Its bound relative to `R` is at most `E_model + E_mesh` under the certified correspondence. Display tolerance may be much larger than 1 nm and must not become model error or evidence for a fine hand-off. This decision does not change the CLI's fields, formats, `exact` flag, legacy mark order or `max_deviation` semantics. Exposing whole-result certificates or stronger export reports through that interface requires a separately reviewed versioned change; no caller may reinterpret the old fields to claim this policy is implemented.

### Field and two-kernel comparison

The field evaluator must retain the same exact recipe revision and mathematical definitions. Its certificate describes its represented regularised solid, independently of B-rep results. A field formula evaluated with enclosing arithmetic can define exact geometry with `E_field = 0`; nonzero evaluation width is query uncertainty, not geometry displacement. Prove that claim for the formula and its domain. A sweep or boolean residual is never a true distance merely because its units are millimetres. The field must supply geometric metric bounds and surviving-trim proofs separately. If it changes geometry through a bounded choice, it owes the same conversion, transfer and publication proofs.

For B-rep and field certificates against the same exact recipe, the triangle inequality gives `D(B, F) <= E_brep + E_field`. A fine hand-off requires this sum strictly below `tau`, or a tighter direct whole-result proof. Two separately passing `0.7 nm` certificates do not establish a sub-nanometre agreement. When comparing the exported file to the field, include `E_export` as well. Exact identities can reduce the total only with proof.

`catallax-compare` must name the target of each comparison. For a query at a mesh sample with unsigned distance enclosure `[l,u]` to the represented field boundary, widening by `E_mesh` gives a local bound against the represented B-rep boundary. To infer distance between the ideal boundaries, also widen by both model certificates:

```text
e = E_mesh + E_brep + E_field
local interval = [max(0, l-e), u+e]
```

Do not widen twice when an input query already encloses distance to the ideal target. Keep an unbounded upper bound unbounded. For a certified whole-domain distance interval `[L,U]` and total threshold `t`, `L > t` proves a discrepancy, `U < t` meets the strict fine threshold, and all other cases remain unresolved. Display comparison keeps its own display threshold and sampled scope. A local upper bound or denser sampling cannot establish whole-result agreement. Fine checks require exhaustive certified coverage, including trimmed faces, missed components, orientation and exact empty results. Origin-set agreement remains a separate obligation; numerical proximity does not prove identity or outward direction.

### Shared algebraic definitions and deterministic work

Implement certificate arithmetic in piece 1 using exact rationals and outward enclosures of shared `Real` expressions. Algebraic dimensions, `RealPoly` coefficients, certified roots and exact poses from the preceding architecture decision supply the definitions and conditioning proofs. Do not fork those representations or rationalise an algebraic frame to simplify the ledger. Keep branch selectors, finite domains and proved positive lower bounds through analytic reduction. A bound on an implicit polynomial requires the same coefficient and domain interpretation in both kernels.

The verification CLI owns allocation, charging, defaults and cancellation semantics for deterministic work budgets. This policy adds proof obligations, not a second budget system. The work owner must account for conversion, transfer and certificate verification as well as the original decision. Exhaustion never enlarges `tau`; it yields a bounded answer only if the required certificate is already established. An interrupted computation cannot publish a partial certificate. No wall-clock observation may select a different geometric branch or certificate. Caches and parallel evaluation must preserve the work owner's deterministic outcomes, including refusal and ordered evidence. Tightening budgets can turn an exact answer into a certified bounded one or a refusal; it cannot justify a wrong answer.

### Observed baseline and evidence

The following are source observations at `b66c5fd73b38eacbe7d2be54de92b8264fa23eac`, not claims that the new policy is enforced. The named tests were inspected for their assertions; this decision does not report a fresh execution of those tests.

| Observed behaviour | Source and existing regression evidence |
| --- | --- |
| `Decider` retains the magnitude of its best enclosure on an exhausted defined decision. `NotExact` has a rational bound and reason but no geometric domain or units | [decider.rs](../crates/catallax-numbers/src/decider.rs), exhausted-limit return; [outcome.rs](../crates/catallax-numbers/src/outcome.rs), `NotExact`; [decisions.rs](../crates/catallax-numbers/tests/decisions.rs) |
| `MarkLog::extend` appends; `worst` takes the maximum per kind. `Answer::is_exact` tests log emptiness | [marks.rs](../crates/catallax-numbers/src/marks.rs); [outcome.rs](../crates/catallax-numbers/src/outcome.rs); [records.rs](../crates/catallax-numbers/tests/records.rs), `marks_worst_by_kind_order_and_answer_mapping` |
| `GeomStore::side_of` converts implicit residuals with `Surface::distance_bound`. Curve incidence can bound distance to carriers without bounding distance to their intersection; normal agreement retains a value bound | [store.rs](../crates/catallax-shapes/src/store.rs); [surface.rs](../crates/catallax-shapes/src/surface.rs); [uncertainty.rs](../crates/catallax-shapes/tests/uncertainty.rs), `intersection_bound_is_to_both_carriers_not_the_curve` and `normals_keep_value_bound_and_decisions_do_not_record_facts` |
| Curve/surface hits can mark a retained unresolved component with `Distance`, an unresolved contact as tangent with `Value`, and merged unordered hits with their point distance | [hits.rs](../crates/catallax-intersect/src/hits.rs), `keep_branch` and candidate contact/order handling; [query.rs](../crates/catallax-intersect/src/query.rs), dimension/point/line conversion helpers. These sites need new whole-effect regressions |
| `Body` stores a mark log. `boolean` appends operation marks then both operand logs; empty paths retain operand marks. The propagation regression injects bounds of hundredths of a millimetre and asserts ordering, not a 1 nm cap | [entity.rs](../crates/catallax-topology/src/entity.rs), `Body`; [lib.rs](../crates/catallax-boolean/src/lib.rs), `boolean`; [sewing.rs](../crates/catallax-boolean/src/sewing.rs), `sew`; [contract.rs](../crates/catallax-boolean/tests/contract.rs), `marks_follow_operation_then_first_then_second_even_with_empty` |
| Measurement accuracy marks contain property enclosure widths as `Value`; `mass_properties` combines its measurement marks without a whole-model displacement calculation | [lib.rs](../crates/catallax-measure/src/lib.rs), `Accuracy::mark` and `mass_properties`; [tests.rs](../crates/catallax-measure/src/tests.rs), `limited_precision_marks_property_widths_in_order`; [measurement.rs](../crates/catallax-measure/tests/measurement.rs), `body_marks_stay_on_body_and_do_not_inflate_measurement` |
| Export constructs an exact answer then adds a distance mark only when `max_deviation` exceeds requested curve tolerance. Its public report has no composed recipe/model bound, and input body marks do not become export marks | [export.rs](../crates/catallax-step/src/export.rs), `Exporter::build`; [lib.rs](../crates/catallax-step/src/lib.rs), `StepReport`; [tests.rs](../crates/catallax-step/src/tests.rs), `rounding_deviation_and_scale` and `branch_spline_has_certified_bounds_at_default_tight_and_small_budgets`; [contract.rs](../crates/catallax-step/tests/contract.rs), `body_marks_are_not_export_marks` |
| The recipe build retains an immutable definition. The CLI separately collects body, measurement and export marks, adds marks when decimal rendering misses measurement/export accuracy, and sets manifest exactness by log emptiness | [evaluate.rs](../crates/catallax-features/src/evaluate.rs), `Built`; [lib.rs](../crates/catallax-cli/src/lib.rs), `evaluate`, `mark_decimal_accuracy` and `mark_entry` |
| Field queries distinguish true distance, residual and unsigned geometric bounds. Comparison widens forward query distances by mesh error and records inherited mesh marks as unresolved | [lib.rs](../crates/catallax-field/src/lib.rs), `Value`, `DistanceBound`, `Query`; [check.rs](../crates/catallax-compare/src/check.rs); [tests.rs](../crates/catallax-compare/src/tests.rs), `unbounded_boolean_upper_distance_never_agrees` and `removed_coplanar_boolean_zeros_stay_unresolved` |

These observations establish why mark propagation alone is insufficient. They do not imply that every `NotExact` is accepted: observed boolean ordering in [graph.rs](../crates/catallax-boolean/src/graph.rs), `sort_points`, refuses ambiguous decisions, and [arc.rs](../crates/catallax-intersect/src/arc.rs) also converts unresolved decisions into refusals.

### Migration and acceptance order

Land interfaces once on the shared algebraic foundation, then migrate producers before accepting their certificates in consumers. Until an operation has migrated, exact paths may remain available with their construction proof; uncertified marked paths must fail closed. No empty-log heuristic may certify an unaudited producer.

1. **Piece 1:** add certificate and typed evidence definitions, rational cap checking and the two refusal codes. Separate diagnostic marks from proof status. Test strict rational comparisons immediately below, at and above `tau`, unit rejection, invalid/missing evidence, deterministic ordering and shared proof identities. Coordinate work accounting with the verification CLI owner without changing that owner's policy.
2. **Piece 2:** implement scoped residual conversions for all shared analytic definitions, algebraic poses and finite domains. Prove lower bounds and branch correspondence. Test scaled implicit equations, small radii, huge lever arms, normals near parallel, singular cone/torus locations and translated algebraic frames. Retain the carrier-versus-curve regression and require an additional curve-distance proof before promoting its result.
3. **Piece 3:** bind certificates to valid bodies and immutable recipe targets, including exact primitives and empty bodies. Enforce topology validity independently of metric closeness. Test forged marks, legacy marked operands, removed components, reversed orientation and multiple touching lumps. Failed construction must not publish a body as successful.
4. **Piece 4:** certify meeting loci, contact classification, finite trims and branch selection. Exercise the observed `hits.rs` mark sites with exhausted budgets. Test near-tangent sphere/plane contacts, nearly parallel planes, clustered hits with unresolved order, open-conic branch changes and narrow finite traces. A tiny polynomial residual with a large meeting displacement must refuse.
5. **Piece 5 and recipe evaluation:** compose through all five stages, each operand and every feature node. Add independent exact-reference checks for three-step accumulation below/at/above the cap, repeated operands, shared subrecipes, chained subtraction, almost coincident cuts, thin cavities and a component that vanishes near tangency. Check operation and body certificates agree. Exact fast paths, empty paths and maximal merging must pass the same publication gate.
6. **Piece 6:** enclose exact-target volume, first moments, centroid and area using model evidence plus numerical integration bounds. Test an analytic offset box with known volume change, thin bodies whose volume interval includes zero, large coordinate offsets and close corrugated patches whose areas differ greatly. Assert intervals enclose independent exact values. A model certificate alone must not license an area or centroid claim.
7. **Piece 7:** compose export with model error and validate the emitted representation. Test a model using `0.8 nm` plus export rounding/fitting using `0.3 nm`, strict equality at 1 nm, individually safe entities with unsafe consecutive transformations, and decimal unit conversion at extreme scales. Read exported values independently and compare certified whole domains, not just writer-side samples. Ensure refusal publishes no usable partial file.

Migrate field, mesh and comparison alongside their dependencies, then require end-to-end recipes with exact, bounded and refused outcomes. Two individually sub-nanometre kernel bounds whose sum exceeds the cap must remain unresolved for hand-off; a certified tighter comparison may settle them. Display agreement must never stand in for the fine proof. Include a component smaller than the sampling spacing and an unsampled slit, plus missing reverse coverage, removed residual zeros and an unbounded distance upper.

Run the same fixtures across thread counts, cold and warm caches, repeat processes and supported machines under identical deterministic allowances. Compare geometry, bounds, proof identities, refusal codes and evidence order, excluding elapsed time. Use analytic/exact references and independently implemented interval checks for correctness; fuzzed samples supplement rather than replace the proof contracts. Threshold tests must use exact rational nanometres, not binary64 equality. Keep adversarial inputs that refuse as unfinished coverage regressions until repaired.

The reassessment point is a complete chain through pieces 1 to 7 and both kernels that demonstrates exact agreement, certified cumulative error below 1 nm, and refusal at or above the cap. Before declaring this policy enforced, inventory every mark producer and consumer, including passthrough and empty-result paths, and show that none bypasses that chain. This ruling authorises implementation of the certificate contract; it does not declare that milestone complete.

## The seven pieces

| Piece | Crate | Job | Depends on |
| --- | --- | --- | --- |
| 1 | `catallax-numbers` | Exact numbers, enclosures, the decision interface, time budgets, construction facts, marks | nothing |
| 2 | `catallax-shapes` | Exact geometric definitions: points, poses, surfaces, curves, vertices by definition; the write-once geometry store | 1 |
| 3 | `catallax-topology` | Faces, edges, loops, shells, bodies over curved geometry; primitive solids; topology checks | 1, 2 |
| 4 | `catallax-intersect` | Exact intersections of plane, cylinder and sphere pairs; curve and surface meetings; closed forms, membership, ordering and tangents along curves; certified curve traces | 1, 2 |
| 5 | `catallax-boolean` | The one boolean engine: intersect, imprint, classify, keep, sew; point containment | 1 to 4 |
| 6 | `catallax-measure` | Volume, area and centre of mass with proven bounds | 1 to 4 |
| 7 | `catallax-step` | STEP AP242 export with a proven deviation report | 1 to 4 |

A front-end crate, `catallax-cli`, builds both `catallax` and `exact-brep` from the same main. It reads a question file and writes a STEP file and manifest. The legacy executable, format strings and `kernel.name` remain unchanged under [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md). An independent tester built without seeing the kernel's code uses that contract.

Build order: pieces 1, 2 and 3 in that order; then 4 and 5; with 6 and 7 alongside them.
Piece 4 needs only 1 and 2, so its work can start as soon as piece 2's interface is stable.

`catallax-recipe` is not a piece of the B-rep half. It depends only on piece 1 for exact numbers. Its public types implement serde's `Serialize` and `Deserialize`, and no serde type appears in its function signatures. Its errors are recipe errors with codes `recipe.<reason>`, not refusals. `catallax-field` uses the shared sketch validator from `catallax-features`, but never calls a B-rep builder, body containment or export.

`catallax-features` evaluates B-rep recipes through the existing operations; field recipe evaluation compiles stock solids, profile sweeps and booleans independently. The two kernels do not convert results after each operation. Conversions happen only when a feature or output needs them.

`catallax-mesh` is a display-only consumer of pieces 1 to 4. It samples each
B-rep edge once, triangulates periodic analytic face charts, and returns shared
positions, face-tagged triangles and proven chord bounds. Its OBJ writer groups
triangles by input face ID. It does not change the body, STEP export or the
question/manifest exchange. Display accuracy is not the PRD's 1 nm agreement
guarantee.

`catallax-compare` consumes both recipe evaluators and the display mesh. It uses
field distance certificates, independent rational boundary samples and exact
trimmed-face projections. It never calls B-rep body containment for comparison.
Its revision controller atomically publishes only complete, current reports.
Regional agreement is limited to the requested box; sampled display agreement
does not certify whole-boundary Hausdorff distance or 1 nm accuracy.

Dependencies point one way only. A piece never depends on a piece below it in the table, except as listed.
Crates outside this list (third-party) may be used inside a piece, but no third-party type appears in any public signature.

## Crate layout

```
Cargo.toml                 workspace, shared lints, no dependencies of its own
crates/catallax-numbers/         piece 1
crates/catallax-shapes/          piece 2
crates/catallax-topology/        piece 3
crates/catallax-intersect/       piece 4
crates/catallax-boolean/         piece 5
crates/catallax-measure/         piece 6
crates/catallax-step/            piece 7
crates/catallax-mesh/            display meshes and OBJ export, beside the pieces
crates/catallax-recipe/          recipe language data model, beside the pieces
crates/catallax-field/           independent analytic fields and certified named queries
crates/catallax-features/        named B-rep recipe builds and sketch-based features
crates/catallax-compare/         two-way sampled display comparison and revision cancellation
crates/catallax-cli/             front end, binaries `catallax` and `exact-brep`
docs/                      this file, DESIGN.md, JUDGE-INTERFACE.md, RECIPE.md, INTERFACE-CHANGES.md, pieces/
```

Every crate forbids `unsafe` code and requires documentation on every public item.
Each crate's `src/lib.rs` re-exports its public types at the crate root, so callers write `catallax_numbers::Real`, not a module path.
Module structure inside a crate is the creator's choice. The re-exported names and signatures are the interface.

## Data flow of one answer

A question names an operation, its operands by exact parameters, and their poses (see `JUDGE-INTERFACE.md`).

1. **Read.** `catallax-cli` parses the question. Every number is parsed as an exact `Rational` (`"1.1"` becomes 11/10). Poses become exact rational rotations and translations.
2. **Build operands.** Piece 3's implemented primitive builders `block`, `cylinder`, `sphere`, `cone` and `torus` add exact surfaces and curves to `GeomStore`, record construction facts and add topology to `Model`. Cone and torus builders do not imply boolean coverage for those surfaces.
3. **Operate.** For each step, piece 5's `boolean` runs its five stages on two bodies of the same `Model`:
   - *intersect* asks piece 4 for every face-pair surface intersection and every edge-face meeting; piece 4 adds intersection curves and vertices by definition to the `GeomStore`;
   - *imprint* splits faces and edges of both bodies along those curves and points, adding new topology;
   - *classify* labels every face piece inside, outside or on the other body, using piece 4's ordering predicates and piece 5's own containment test;
   - *keep* selects the face pieces the operation keeps;
   - *sew* joins them into a new body.

   Every decision on the way is made by piece 1's `Decider`, which consults the store's construction facts first (through piece 2), then the fast filter, then more precision, then exact methods. The result is a new `BodyId`; the operands stay valid.

   When both arguments of a step resolve to the same `BodyId`, `catallax-cli` applies the set rule without calling `boolean`, which forbids identical body handles. Union and intersection return that body; subtraction creates a new `Empty` body. Each result carries the input body's marks once, and the step makes no decision.
4. **Measure.** Piece 6 computes volume, area and centre of mass of the result as proven enclosures, using piece 2's exact surface parametrisations and piece 4's certified curve traces.
5. **Export.** Piece 7 writes the body as STEP AP242. Analytic geometry is written as analytic entities. Curves with no exact STEP form are approximated, and the largest distance between the exported and the exact geometry is proven and reported.
6. **Report.** `catallax-cli` writes the manifest: status, measured properties with bounds, every mark with its reason and bound, any refusal with its piece and stage, the export deviation, and the measured work and time.

## The decision outcome vocabulary

Defined in piece 1 and used by every piece. Two levels exist.

### Decision level: `Outcome<T>`

Every single decision returns `catallax_numbers::Outcome<T>`:

| Variant | Meaning | What the caller does |
| --- | --- | --- |
| `Proven(T)` | The answer `T` is certified; the quantity behind it is proven non-zero. | Uses `T`. |
| `ProvenCoincident(Coincidence)` | The quantity is proven exactly zero. `Coincidence` says how: exact arithmetic, an algebraic proof, a separation bound, or a recorded construction fact. | Takes the degenerate branch (touching, tangent, coincident, on the boundary). |
| `NotExact(NotExact)` | The budget ran out. The quantity is proven to lie in `[-bound, bound]`; its sign is unknown. | Takes the degenerate branch as if coincident, and records a `Mark` built from the `NotExact` on the result. |
| `Refused(Refusal)` | The decision cannot be answered (for example a square root of a proven-negative number, or a case outside the piece's coverage). | Stops the operation and passes the refusal up. |

`T` is the strict answer: `StrictSign` (`Negative`, `Positive`), `StrictOrder` (`Less`, `Greater`), `Side` (`Front`, `Back`), `Containment` (`Inside`, `Outside`) and so on. The degenerate answer is never a value of `T`; it is always the `ProvenCoincident` variant. This keeps every caller honest about the degenerate case: the compiler makes it handle it.

Observed inherited behaviour uses `NotExact` as a bounded, visible tolerance at callers that accept it. Under the composable geometric error decision above, taking that branch requires the scoped geometric and whole-result proofs. No caller may instead pick a side from the cached approximation.

The table describes observed inherited behaviour, not completion of the new PRD. The composable geometric error decision above defines the required conversion, composition and strict 1 nm cap. Enforcement remains implementation work; this publication changes no numerical code or outcome API.

Only impossible geometry is an acceptable final refusal. Unsupported pairs, ambiguous classifications, exhausted budgets, deadline refusals and missed measurement accuracy are unfinished implementation work. Preserve their honest stable codes until the underlying coverage is repaired.

### Operation level: `Answered<T>`

Every operation (an intersection, an imprint, a boolean, a measurement, an export) returns `catallax_numbers::Answered<T>`, which is `Result<Answer<T>, Refusal>`.
`Answer<T>` holds the value and a `MarkLog`: every `NotExact` met on the way, each turned into a `Mark` with its reason code, its bound and the kind of bound (a distance, an angle, or the raw quantity).
Observed inherited `Answer::is_exact` treats an empty mark log as exact. Observed marked answers retain local bounds, and `Body` stores the marks of the operation that produced it, merged with the marks of its operands. Under the new policy, exactness and whole-result accuracy require a certificate independently of log emptiness; a marked measurement can still describe exact model geometry.
The boolean refuses with `BooleanRefusal`, which adds the `Stage` that broke, so failures can be filed by the first stage that broke.

### Refusal and mark codes

`Refusal::code` and `Mark::code` are stable strings of the form `piece.reason`, for example `intersect.unsupported-pair` or `numbers.budget-exhausted`.
The piece prefixes are `numbers`, `shapes`, `topology`, `intersect`, `boolean`, `measure`, `step` and `cli`.
Each piece description lists its codes. Adding a code is an interface change and is logged.
`detail` is free text for people and is not part of the contract.

## The time budget

Each decision runs under a budget set by `catallax_numbers::BudgetPolicy`:

- `max_precision_bits`: the most bits of working precision the decision may use while refining;
- `max_work`: the most work units the decision may spend in total.

A work unit is defined by piece 1 (see its description) as a deterministic count of limb operations, so the same decision spends the same work on every machine.
A `Decider` holds the policy and creates a fresh budget for every decision. It also keeps `DecisionStats`: how many decisions ended on each proof path, how many were marked or refused, the largest precision used and the total work.
The policy default is set by the architect, recorded in `INTERFACE-CHANGES.md` when it changes, and passed by the front end. Pieces never invent their own budgets.

At the level of a whole operation, the `Decider` may carry a wall-clock deadline. A decision that starts after the deadline, or is still in progress when it passes, returns `Refused` with code `numbers.deadline`. This is the only way wall-clock time affects the kernel, and it can only produce a refusal, never an answer.

The `Decider` checks the deadline during budgeted work and before returning a decision. `check_deadline()` lets callers abandon long arithmetic, subdivision and integration loops without starting a numeric decision or changing its statistics. Trace enclosures use the caller's `Decider` and propagate deadline refusals. Fact-only piece operations may still answer without a deadline check; the CLI checks the whole answer before publishing it and removes any partial STEP file on a deadline refusal.

Decision of 2026-10-22: publish `catallax-verify` and `catallax-request/1` as a separate deterministic exchange. Each request receives 100,000,000 work units by default, adjustable as explicit input. `Decider` and `Work` checkpoints cost one unit and accepted arithmetic reservations retain their existing work cost. All deciders in one synchronous request share the allowance. Exhaustion latches `numbers.operation-budget` and prevents publication of partial output. Worker count changes request scheduling only. Wall time is reported separately and never selects an answer, refusal or mark. External cancellation is an incomplete run, not a kernel outcome. The legacy commands retain their 50-second behavior for compatibility. See [VERIFY-INTERFACE.md](VERIFY-INTERFACE.md).

Each request runs on a fresh thread inside the persistent process. The worker-count proof found that thread-local intersection plans carried expression sharing across requests and changed work totals. Thread teardown isolates that cache history without changing geometry or the legacy commands. Numeric execution remains synchronous inside the request. The measured wall time includes thread startup and teardown.

The default allowance is a finite safety policy, not a latency promise or a conversion from 50 seconds. It allows substantial exact work while bounding repeated decisions that individually fit their budgets. The allowance and accounting rule are versioned inputs to verification; changing the default requires a contract review. A time limit, CPU-time limit or load-calibrated allowance would still let machine behavior choose outcomes, so none is used.

Interactive performance targets are 100 ms for stock solids, 500 ms for booleans, 1,000 ms for recipe feature builds and 2,000 ms for the complete background two-kernel comparison. Stock creation should feel immediate. Booleans and feature edits should update within a short interaction; the more expensive sampled comparison runs in the background. These targets include requested measurement and export work, making costly export visible rather than hiding it. They are engineering targets, not correctness gates. Slow cases are performance bugs even when their proofs are correct. Measure on a pinned core with warmup and matched AB/BA pairs; report noisy shared-machine results as such.

Manifest/3 adds distinct `elliptic-cylinder`, `parabola` and `hyperbola` kind tags. Recipe and comparison adapters consume the same shared recipe definitions. Existing tags keep their meanings. This closes the separate exchange decision required by the exact geometry ruling; implementing missing geometry remains other work.

## Recorded construction facts

Piece 1 defines the vocabulary and the store: `FactKind` (coincident, coaxial, concentric, parallel, tangent, incident), `Fact<K>`, `FactOrigin` (construction or proof) and `FactBook<K>`, generic over the key type.
Piece 2 owns the geometry keys (`GeomKey`: a surface, curve or point id) and holds one `FactBook<GeomKey>` inside its `GeomStore`.

- Builders record facts when the construction guarantees them. Piece 3's cylinder builder records that each cap circle is incident on its cap plane and on the cylinder. Piece 4 records that every intersection curve it creates is incident on both of its surfaces, and that every vertex it creates is incident on the geometry it was solved from.
- A decision that proves an exact zero by numbers may record that zero as a fact with origin `Proof`, so later decisions reuse it.
- Piece 2's geometric decisions that name stored objects (`side`, `point_on_curve`, `same_point`, `same_surface`) look up facts first. `side_of` takes a bare `Point3` and has no key to look up. A matching fact yields `ProvenCoincident` with `ZeroProof::ConstructionFact` without any arithmetic.
- Equivalence kinds (coincident, coaxial, concentric, parallel) are answered transitively by the fact book. Incidence chains (a point on a curve that lies on a surface) are answered transitively too.

## Write-once storage, handles, coordinates as a cache

Piece 2 defines `Handle<T>` (a typed index) and `Arena<T>` (an append-only store).
Once an entity is added it never changes. There is no `get_mut`.
Operations create new entities and return new handles; the operands of a boolean stay valid and unchanged afterwards, so a chain of booleans can be replayed or compared step by step.

- `GeomStore` (piece 2) holds surfaces, curves and point definitions in three arenas, plus the fact book.
- `Model` (piece 3) owns a `GeomStore` and holds vertices, edges, coedges, loops, faces, shells, lumps and bodies in their own arenas.

Coordinates are a cache.
A vertex stores a `PointId`, and the point is a definition (`PointDef::Explicit` for an exact rational point, `PointDef::Meet` for "where these meet").
The exact position is available as `Point3` with `Real` coordinates, from which enclosures of any precision can be computed.
Cached `f64` coordinates and cached enclosures may be kept next to an entity for speed, but they are derived data: dropping them changes no answer.
Handles are only meaningful within the store that issued them. Mixing stores is a programming error.

Storage is single-writer for now: operations take `&mut Model`. Read-only queries take `&Model`. Multithreaded execution comes later and must keep the determinism rule.

## Recipe names and two-kernel agreement

Future persistent face and edge names identify the creating operation and its inputs, never the traversal order or an arena index. Store-local `FaceId` and `EdgeId` handles remain useful during evaluation but are not rebuild-stable recipe references. [RECIPE.md](RECIPE.md) section 3 defines the splitting, merging and feature-edit naming rules; no kernel reports these names yet.

Both evaluators use the same operation definitions and agree exactly where mathematics permits. The only necessary approximations are those in PRD section 3.1, bounded by 1 nm. Definition-exact surfaces that STEP cannot express are approximated only on export.

Planned background checks compare the field to the B-rep display mesh near an edit. This coarse check cannot certify 1 nm. Development checks and finer whole-part checks at output hand-offs need separate proofs and measured error bounds. An agreement failure fails the edited feature; the application must never show a half-built body as a successful result.

Measure rebuild, query and comparison latency per edit. Safe Rust threads may rebuild and check in the background once determinism is preserved. The current `&mut Model` evaluation path remains single-writer.

## Faces and edges carry curved surfaces from day one

Piece 3's types never assume planarity:

- `Face` holds any `SurfaceId`, a sense relative to the surface normal, and any number of loops, including none (a whole sphere) and several of equal standing (a cylinder band between two circles). There is no "outer loop" field.
- `Loop` is either a cycle of coedges or a single vertex (the apex of a cone).
- `Edge` holds any `CurveId`, a sense relative to the curve's natural direction, and optional start and end vertices. An edge with neither is a ring edge covering a whole closed curve, such as the rim of a cylinder.
- `Vertex` holds any `PointId`, so a vertex can be the meeting point of three curved surfaces.
- Coedges are separate from edges, and an edge may have more than two coedges, so imprinting and non-manifold results fit in the same structure.

Seams are a property of a surface's parametrisation, not of the topology. A cylinder face needs no seam edge. Poses "along the seam" are therefore a test of parametrisation-dependent code (measurement and export), not of topology.

## Boundary types by crate

The full signatures are in the piece descriptions and the crate skeletons. This table lists what crosses each boundary.

| From | To | Types |
| --- | --- | --- |
| 1 | all | `Integer`, `Rational`, `Sign`, `Interval`, `Ball`, `UPoly`, `MPoly`, `Algebraic`, `RealRoot`, `Real`, `PolySystem`, `RootBox`, `IsolatedRoot`, `BoxRoots`, `Outcome`, `StrictSign`, `StrictOrder`, `Coincidence`, `ZeroProof`, `NotExact`, `Refusal`, `Answer`, `Answered`, `BudgetPolicy`, `Decider`, `DecisionStats`, `FactKind`, `Fact`, `FactOrigin`, `FactBook`, `Mark`, `BoundKind`, `MarkLog` |
| 2 | 3, 4, 5, 6, 7 | `Handle`, `Arena`, `RPoint3`, `RVector3`, `Point3`, `Vector3`, `Rotation`, `Pose`, `Surface` and its kinds, `SurfaceKind`, `ParamKind`, `Curve` and its kinds, `CurveKind`, `PointDef`, `Meet`, `MeetOf`, `Position`, `GeomStore`, `GeomKey`, `SurfaceId`, `CurveId`, `PointId`, `Side`, `Apart` |
| 3 | 5, 6, 7, cli | `Model`, `Vertex`, `Edge`, `Coedge`, `Loop`, `Face`, `Shell`, `Lump`, `Body`, `BodyKind`, `Sense`, their ids, `EntityCounts`, `TopologyFault`, builders `block`, `cylinder`, `sphere`, `cone`, `torus`, `empty` |
| 4 | 5, 6, 7 | `SurfacePair`, `Component`, `Contact`, `CurveSurface`, `CurveHit`, `Cyclic`, `CurveForm`, `CurveTrace`, `TracePiece`, `intersect_surfaces`, `intersect_curve_surface`, `curve_form`, `curve_is_closed`, `on_curve`, `same_curve`, `order_on_curve`, `cyclic_order`, `tangent_at`, `arc_point`, `trace`, `supports` |
| 5 | cli | `BoolOp`, `Stage`, `StageRecord`, `BooleanAnswer`, `BooleanRefusal`, `StageReport`, `SectionEdge`, `IntersectionGraph`, `Imprinted`, `FaceClass`, `Classification`, `Selection`, `Containment`, `boolean`, `intersect`, `imprint`, `classify`, `keep`, `sew`, `classify_point`, `face_point` |
| 6 | cli | `Enclosure`, `Accuracy`, `MassProperties`, `mass_properties`, `volume`, `area`, `face_area` |
| 7 | cli | `StepOptions`, `StepReport`, `write_step` |

There are no public traits at the piece boundaries in the skeleton. Pieces call each other through concrete types and functions.
This is deliberate: a trait seam would invite a second implementation of a piece to diverge silently in behaviour.
The comparison step of the blind build adapts an existing equivalent to these concrete signatures instead.
The one exception planned for later is meshing for display, which is plumbing and may sit behind an adapter.

## Determinism and order

Every collection that a result depends on is ordered: `Vec` in a defined order, or `BTreeMap`/`BTreeSet`.
Hash maps may be used for lookup only, never iterated to build a result.
Where a piece returns several items (intersection components, hits along a curve, faces of a body), its description defines the order.

## Errors that are not refusals

Passing a handle from another store, an index out of range, or a malformed polynomial system is a programming error and may panic.
Anything that can come from a question file (zero-size blocks, zero quaternions, unsupported shapes, unsupported pairs) must be a `Refusal`, never a panic.

## Changing this architecture

Workers propose interface changes to the architect. The architect decides, updates the affected description and crate skeleton, appends an entry to `INTERFACE-CHANGES.md` naming the date and the affected pieces, and announces it to every live worker.
Nobody else edits a public signature.
