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

Amended by the owner on 2026-10-23. Certify the whole result as exact except for its listed approximated or imported faces and recorded bounded decisions. Each exception carries a proved local tolerance and its applicable limit. This replaces the original universal requirement that the represented solid and every hand-off be strictly below 1 nm from the recipe. The 1 nm modelling tolerance remains the limit for sameness and decisions. A fitted face within its local cap does not also need a sub-nanometre whole-solid certificate.

A mark records why a decision stopped; it is not a geometric proof. This is the implementation contract for [PLAN.md](PLAN.md)'s error-policy item and [PRD section 3.1](PRD.md#31-agreement). It supersedes the inherited rule below that automatically takes every `NotExact` as coincidence. This documentation decision changes no code, work-budget policy or existing exchange contract. An unchanged exchange that promises 1 nm must retain that promise or refuse; exposing the amended certificate needs a versioned contract.

| Decision | Reason | Alternative and cost |
| --- | --- | --- |
| Certify a complete partition of exact regions and explicit local exceptions against the recipe | Every face, trim and material region remains accounted for | Certify only individual predicates or sample points. Small residuals can select a distant branch or remove material |
| Keep 1 nm for sameness and decisions; accept fitted/imported faces under their recorded local tolerances and caps | The permitted approximation belongs to a named piece with an exact definition | Require every fit and output to meet a universal 1 nm bound. This defeats the approved local-tolerance policy |
| Retain typed evidence and a shared dependency graph | A consumer can verify scope, units, conditioning and inherited effects | Sum or take the maximum of the mark log. Neither operation supplies missing geometric proofs |
| Require a transfer proof or direct remeasurement for each operation | Boolean boundaries and topology can change discontinuously near contact | Assume booleans are non-expansive because field expressions use `min` and `max`. Those expressions need not be true distance |
| Keep property, display and export approximation bounds separate | Area, volume and output chords have different meanings | Treat every bound as a modelling tolerance. This accepts dimensionally invalid evidence and hides output error |

### What the whole-result certificate proves

Let `R` be the regularised solid defined by the immutable exact recipe prefix, including its branch and domain rules. Let `B` be the valid solid represented by the returned body. Neither is a mesh or a numerical sample. The certificate covers the whole result, but its acceptance limits are local:

```text
tau = 1/1_000_000 mm
whole-result certificate = validity + complete correspondence + material orientation
                         + exact regions + explicit local exceptions + decision records
exception i = source, exact definition, represented domain, U_i, applicable cap, proof
```

For each exception, `U_i` is a nonnegative exact rational upper bound on deviation over the complete corresponding trimmed domain. Cover both directions, including its edges and vertices; prove the branch and outward material orientation. The certificate proves exact agreement on the remaining geometry and confines any material disagreement to the certified uncertainty regions. Boundary Hausdorff distance alone would miss an inside/outside reversal. A list of face tolerances without these coverage and material proofs is not a whole-result certificate.

| Exception source | Acceptance rule |
| --- | --- |
| Kernel-made shape approximation | Target `U_i <= 100 nm`; hard cap `U_i <= 1 µm` |
| Imported geometry | `U_i <= 10 µm`, or the explicit looser limit accepted and recorded for that file |
| Bounded decision only | Prove the substituted geometric gap strictly below `tau` and its effect over the affected domain. A decision-only geometric exception must remain strictly below `tau`; it cannot borrow a fit/import allowance |

Every local exception also respects the tenth-of-the-shortest-neighbouring-edge and tenth-of-the-narrowest-neighbouring-face limits in PLAN section 3.4. Record proved size bounds. The kernel-made target is an engineering target; exceeding it while remaining within every hard cap is reported, not falsely called a cap failure. Imported source uncertainty and deviation from the imported file's declared geometry are distinct quantities with distinct targets. Retain both when provided; do not invent an exact intended design behind a file.

An approximated face's allowance does not enlarge the gap allowed for a later sameness decision on it. Record that decision separately and remeasure the resulting face against its retained exact definition. Do not relabel a budget-exhausted decision as a shape with no exact carrier to gain the larger cap. A face whose support is exact but whose trim moved still needs an exception for that changed domain.

Two empty results agree exactly. A nonempty recipe result cannot be certified as empty by deleting an exception. Prove exact emptiness or refuse. Assemblies need per-body certificates and a matching proof; overlapping bodies are not implicitly replaced by their union. No removed component, hidden cavity or unlisted face may disappear into a local tolerance.

For consumers that need a uniform distance bound, the certificate may additionally prove

```text
D(B,R) = sup over x in Euclidean 3-space of |s_B(x) - s_R(x)| <= E
```

Here `s_S` is true Euclidean signed distance, negative inside. This summary needs a complete oriented correspondence proof; taking `max(U_i)` without that proof is invalid. There is no universal requirement `E < tau` for bodies containing permitted face approximations. The local records remain authoritative even when a uniform summary exists. Proven zero displacement means exact geometry; an empty mark log, zero cached residual or exact-looking file does not. Query enclosures may have width even for exact geometry.

The certificate identifies the recipe revision, represented body, units, target definitions, covered domains, bounds and proof dependencies. Preserve each decision code, residual enclosure, affected entities and branch choice. Use stable definition identities and ordered proof records, never pointer addresses or worker completion order. Dependencies form an acyclic graph. Geometry remains stored by definition; approximations are explicit representations beside it. A mark never authorises a false construction fact or an invalid shell whose edge merely lies near its face.

### From a local mark to millimetres

Every geometry-changing use of `NotExact` needs a conversion and a proof of the effect of taking the degenerate branch on the entire affected domain. Each accepted coincidence substitute must have a proved geometric gap strictly below `tau`. Account for its consequences in the relevant local records and the whole-result cover. Later cancellation cannot excuse an over-limit local choice. A conversion at one witness does not cover a face, trimmed curve or body. Piece 1 cannot infer geometric units from a `Real` expression. The geometric caller proves every denominator, domain restriction and conditioning bound used in the conversion.

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

For the exact operation `F`, represented inputs `B_j`, exact inputs `R_j` and each affected output region `i`, prove a scoped transfer relation or certify the output directly against `F(R_1, ..., R_n)`:

```text
U_out,i <= T_F,i(input local certificates and domains) + delta_F,i
```

`delta_F,i` covers new representation changes. `T_F,i` covers inherited effects on that output region, including which material survives. Each resulting bound must meet its own applicable cap; the output certificate must cover all regions. Neither a default maximum nor a sum of input tolerances proves a transfer for arbitrary booleans. Near contact, a small input displacement can create or destroy a component or change which distant patch survives.

Consecutive representations of the same target compose by the triangle inequality unless a tighter joint proof is available. A Lipschitz transfer needs proved constants over the whole uncertainty domain. Direct remeasurement against the retained exact definition can replace a loose accumulated estimate, with its dependencies retained. A fresh operation does not grant a fresh displacement allowance relative to yesterday's approximation. Independent face allowances are not summed into a universal 1 nm budget.

| Stage | Required proof before its result can be used |
| --- | --- |
| Intersect | Retain every relevant meeting, component selector and finite trim. Convert residuals to bounds for the intended meeting locus. Prove the effect of any omitted or merged meeting, including downstream classification |
| Imprint | Prove that the split arrangement covers each input face, cyclic ordering and incidence are valid, and merged cuts stay within their certified correspondence. Carry intersection evidence forward without counting it twice |
| Classify | Certify material labels outside the recorded uncertainty regions, with oriented evidence inside them. A tangent or on-boundary branch must not select an unbounded or distant face. An unresolved sample cannot label a whole face |
| Keep | Derive the selected material from regularised set semantics for both operands. Account for lost and retained patches, thin walls, cavities and emptiness. A small local gap does not prove a component dispensable |
| Sew, including maximal merging | Prove a valid solid and complete oriented correspondence. List every approximated/imported face and bounded decision, verify its local limit and every exact region, then publish the whole-result certificate |

Intermediate stage objects carry unresolved obligations explicitly and cannot be exposed as certified bodies. Check applicable caps whenever enough evidence exists and again at publication. A stage may retain alternatives for later proof; it may not discard an unresolved alternative and present the result as certified.

Recipe evaluation follows the same rule at every body-producing node. Exact rigid motions preserve local bounds and transport their scopes; a certified uniform scale multiplies lengths and bounds by the absolute scale, then rechecks the dimensional caps. Other features need their own transfer or direct proof. Reusing a body in two operands shares proof identity but does not remove either influence. Deduplicate evidence, never numeric contributions by matching reason strings. Exact identities such as `A subtract A = empty` can discharge inherited displacement for that output while retaining the audit record. Equal-looking bodies need an equality proof first.

There is no blind per-operation reset. For example, a 50 nm fit remeasured at 60 nm against its original definition records 60 nm, not a fresh 10 nm. It can pass its local caps without being refined below 1 nm. A chain of bounded identifications still needs evidence against the shared representative; three individually small gaps do not establish transitive sameness. If composed evidence exceeds a region's limit, refine or reprove it directly; never clamp the number or hide it in another face's allowance.

### Refusals and publication

Reserve these stable codes in piece 1. The implementation entry identifies the producing piece, recipe step, affected region, limit and boolean stage where applicable.

| Code | Trigger |
| --- | --- |
| `numbers.geometric-gap-unproven` | A local conversion, transfer theorem, valid representation, complete scope or whole-result certificate is missing or cannot be proved within the allowance |
| `numbers.geometric-gap-exceeded` | The best certified decision gap is at least `tau`, or a certified local upper bound exceeds its applicable hard cap or requested hand-off limit |

The second code says the available proof does not meet the named limit, not that the true error necessarily exceeds it. Decision acceptance is strict; face hard caps are inclusive. Refine, split/refit or take an exact branch within the assigned work allowance if possible. Otherwise return the specific refusal and no successful body or output for the failed operation. Do not clamp bounds, enlarge caps without the authorised import acceptance, drop evidence or call a bounded coincidence `ProvenCoincident`. Proof or budget failures remain unfinished coverage under PRD section 3.2. A proved self-touching result that violates the solid rules remains an impossible-geometry refusal; a singular chart alone does not establish one.

Every piece rejects legacy marked input without sufficient evidence, preserves certificates through passthrough paths and validates imported or deserialised evidence. A producer cannot self-certify by constructing a public mark. Final publication checks the complete certificate, including changes after the last boolean. Failure to meet a requested output accuracy does not retrospectively invalidate a model that passed its own local caps.

### Measurements, STEP and display

Piece 6 states whether an enclosure describes `B` or the ideal recipe solid `R`. Shared mass properties must enclose `R`. Numerical quadrature and rounding error remain separate from model displacement. Integrate the exact definition directly, or use the certificate's complete uncertainty regions to bound the difference between `B` and `R`. A proved upper volume for the union of their local tubes bounds the volume difference; do not substitute `area*U` near corners or contacts. Bound coordinate magnitudes in those regions to bound first-moment differences. Divide by a volume interval proved positive to enclose the centroid; otherwise preserve `measure.centroid-not-separated`.

A small positional tolerance alone gives no general area bound. Area needs certified correspondence and Jacobian bounds or direct exact-definition integration. Missing evidence refuses that property. A 50 nm face tolerance is compatible with model acceptance; a requested property interval may still require tighter fitting or direct integration. This is a property-accuracy request, not a universal sub-nanometre model gate. Diagnostic intervals that miss requested accuracy must not be reported as successful hand-offs.

Piece 7 certifies the solid described by the emitted surfaces, numbers, curves and trims. Give each exported region an output-displacement bound `V_export,i` relative to its model region. Under the certified correspondence, report its total bound against the exact definition:

```text
W_file,i <= U_model,i + V_export,i
```

Exact model regions have `U_model,i = 0`. Include decimal rounding, fitting, endpoints, edge/face compatibility, subdivision and added seams in the export proof. Preserve the model's local record unchanged and identify any file-only faces that correspond to it. A file must meet the declared export accuracy, whose contract states whether it measures file-to-model displacement or total file-to-definition deviation. For a requested total bound `A_i`, prove `W_file,i <= A_i` or a tighter direct bound. There is no automatic `U_model,i + V_export,i < 1 nm` requirement. An old or explicitly requested strict 1 nm output contract still requires it; use a versioned exchange for the new policy rather than silently weakening an existing promise.

The default output report must declare its accuracy, target and composed bounds; it must not describe a file as within the model face's local tolerance if rounding made that claim false. A request to retain that same tolerance in the file needs a proof within it, tighter output construction or a specific export refusal. Export settings do not enlarge the model's local cap or rewrite its evidence. Check the full file certificate before publication; a refused export leaves no partial artifact presented as usable.

Display meshing similarly retains a separate `V_mesh,i`, including rounding and complete triangle coverage. Its bound against the definition is at most `U_model,i + V_mesh,i` under the proved correspondence. Display accuracy may be coarser than modelling or export accuracy and never substitutes for either certificate. The existing CLI fields, `exact` flag, mark order and `max_deviation` semantics remain unchanged until a separately reviewed versioned extension exposes these records.

### Field and two-kernel comparison

The field retains the same exact recipe revision and definitions. A formula evaluated with enclosing arithmetic can represent exact geometry with zero model displacement; nonzero evaluation width is query uncertainty. Prove equality for the formula and domain. A sweep or boolean residual is not true distance merely because its units are millimetres. Metric and surviving-trim proofs remain separate. Where the field takes over a B-rep result for an unsupported step, it retains that step's local exceptions and reports the loss of independent construction coverage.

For corresponding regions of independently certified B-rep and field results, the triangle inequality gives a comparison allowance `U_brep,i + U_field,i`. An exact field definition contributes zero, so the B-rep face is checked against its recorded local tolerance. Verify each kernel's own bound against the shared definition as well as their mutual agreement; two kernels can agree on the same wrong shape. An output comparison adds its scoped mesh/export bound only when the query actually targets that output. Do not widen twice when a query already encloses distance to the exact definition.

For a complete certified discrepancy enclosure `[L,U]` and the applicable inclusive local allowance `t`, `L > t` proves failure, `U <= t` proves the metric part of acceptance, and otherwise the result is unresolved. A decision's strictly-less-than-1-nm rule retains its strict comparison. Equality of exact regions needs an exact proof, not a positive comparison allowance invented for convenience. An unbounded upper end never establishes closeness.

Fine whole-part hand-offs check complete coverage, trims, branches, orientation, origin sets, empty results and every local limit. They no longer demand a universal 1 nm difference on permitted approximated faces. A 50 nm face with a complete proof within its recorded tolerance can pass; one exceeding that tolerance cannot. A sample comparison remains a display check regardless of sample density. Include missed components, small slits, removed residual zeros and missing reverse coverage in the rejection tests.

### Shared algebraic definitions and deterministic work

Implement certificate arithmetic using exact rationals and outward enclosures of shared `Real` expressions. Algebraic dimensions, `RealPoly` coefficients, certified roots and exact poses supply definitions and conditioning proofs. Do not fork them or rationalise a frame to simplify the records. Retain branch selectors, finite domains and proved positive lower bounds through analytic reduction. An implicit residual requires the same coefficient and domain interpretation in both kernels.

The verification CLI owns allocation, charging, defaults and cancellation. Account for conversion, transfer, local measurement and certificate verification as well as the original decision. Exhaustion enlarges neither `tau` nor any local cap. Interrupted work cannot publish a partial certificate. Wall time cannot choose branches or acceptance bounds. Caches and parallel work must preserve deterministic outcomes, evidence order and logical work allowances. A tighter budget may produce a certified bounded answer or a refusal; it cannot justify an unproved answer.

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

Land interfaces once on the shared algebraic foundation, migrating producers before accepting their certificates. Keep audited exact paths usable; uncertified marked paths fail closed. Each tolerance-core piece lands with the feature it unlocks, starting with the approximated offset face specified below. Certificate definitions alone are not that milestone.

1. **Piece 1:** add typed local evidence, source/limit records, the complete result certificate and refusal codes. Test decisions below/at/above 1 nm; kernel fits at 50 nm, 100 nm and the inclusive 1 µm cap; imports at 10 µm and a recorded override; feature-size caps; forged evidence, units, order and identity. A valid 50 nm fit must not fail because it exceeds 1 nm.
2. **Piece 2:** implement scoped conversions and the selected spline carrier's certified charts, domains and derivatives. Test scaled residuals, small radii, large lever arms, near-parallel normals, singularities and algebraic poses. Retain the carrier-versus-curve regression; it needs a curve correspondence proof.
3. **Piece 3:** bind local records and complete coverage to valid bodies and recipe targets. Test missing exceptions, moved trims on exact supports, wrong orientation, removed components, forged incidence, empty results and stale revisions. Test an exact region beside a valid 50 nm approximated face without granting the exact region 50 nm.
4. **Piece 4:** certify meeting loci, contact, finite trims and branch selection. Test near-tangency, clustered roots, open branches and narrow traces. Even on a fitted face, a 2 nm sameness gap fails the 1 nm decision rule. A tiny residual with an unbounded geometric effect remains unproven.
5. **Piece 5 and recipe evaluation:** compose scoped evidence through every boolean stage and feature node, then remeasure against exact definitions. Test repeated operands, changed trims, a second offset, thin cavities, shrinking feature-size limits and nontransitive identifications. Test a 50 nm fit rebuilt at 60 nm with that actual bound recorded rather than reset to zero or forced below 1 nm.
6. **Piece 6:** enclose target volume, moments, centroid and area using exact integration or local uncertainty evidence plus quadrature bounds. Test thin bodies, large translations and corrugated close surfaces. A valid model can fail a tighter property request without being declared invalid.
7. **Piece 7:** certify export bounds separately and compose them per region. A 50 nm model bound plus 20 nm export error can establish a 70 nm total bound under a proved correspondence; it cannot establish a 50 nm total export request. Test stricter requested output, decimal scale extremes, edge/face compatibility, new seams and independently read whole domains. Existing strict-output exchanges retain their promises.

Migrate field, mesh and comparison alongside their dependencies. Test exact, locally bounded and refused end-to-end recipes, including an exact field checking a 50 nm fitted face and two fitted representations whose mutual allowance is the sum of their separately proved bounds. Neither a uniform summary nor a matching sampled mesh replaces the complete certificate. Run identical fixtures at multiple thread counts, cold/warm caches and repeated processes under the same deterministic allowance; compare evidence and refusals, excluding elapsed time.

Completion requires a feature chain through both kernels, with complete exact/exception coverage, local caps, valid topology, scoped property/output hand-offs and correct refusals. Inventory every producer and consumer, including passthrough and empty paths, before declaring the policy enforced. This ruling authorises implementation; it does not claim that enforcement has landed.

## Faces beyond the analytic kinds: ruling of 2026-10-23

Adopt a certified polynomial or rational B-spline surface beside the immutable exact feature definition, option A below. The spline is the represented geometry of the B-rep face. The definition remains the target for construction, error proofs and later features. Keep existing analytic reductions exact. For the first offset, use a cubic spline cross-section extruded exactly along the cylinder axis, a tensor-product surface of degrees 3 and 1. This also gives general fillets and blends a common candidate representation without making triangles the model.

The owner approved carrier A and the PRD case-1 widening on 2026-10-23, together with the amendment to [the composable error ruling](#composable-geometric-error-decision-of-2026-10-22). The whole-result certificate now records exact regions and every permitted local exception. A fitted face within its recorded tolerance and caps does not also need refinement below a universal 1 nm bound. This docs-only change supplies the representation contract; construction, certificate enforcement and measuring-stick extensions remain implementation work.

The owner decided the geometry-creation tolerance at 1 nm on 22 October, as PLAN sections 3.1 and 15 record. Keep that value for sameness and decisions. Local approximation caps, property accuracy and display/export settings have their own declared scopes; none enlarges the decision gap.

### Why the first offset needs a new carrier

For an ellipse with semi-axes `a > b > 0`, its outward offset by `d > 0` has support function

```text
h_d(theta) = sqrt(a² cos²(theta) + b² sin²(theta)) + d.
```

It is not the ellipse with semi-axes `a+d` and `b+d`. At `a=2`, `b=1`, `d=1` and 45 degrees, that ellipse has squared support `13/2`; the true offset has squared support `7/2 + 2 sqrt(5/2)`, which is larger since `10 > 9`. The extrema along the symmetry axes fix the only possible axis-aligned ellipse, so this also excludes another ellipse with different radii. Extruding the curve does not turn it into any of the six current analytic kinds. This is an exact distinction.

The exact definition is still available. On a regular patch it is `X_d = X + d*n`, with the chosen unit normal, finite domain and branch. For the ellipse, write `X=(a*c,b*s,z)` and `n=(b*c,a*s,0)/sqrt(b²*c²+a²*s²)`, where `c²+s²=1`. Rational stereographic charts for `c,s` give algebraic expressions on bounded chart intervals; no sampled point cloud defines the offset. Near focal points or changing nearest witnesses, this normal formula needs additional domain and envelope proofs. It is not a definition of the entire offset of an arbitrary trimmed solid. [RECIPE section 4.8](RECIPE.md#48-offset) supplies that true-distance definition.

### Carrier comparison

A, B and C all retain the exact operation definition, input revisions, distance or radius, selected branches, surviving domains, material orientation and source origins. None may replace that record with fitted coordinates. The following are design requirements and cost estimates, not observed timings.

| Consumer | A: certified spline beside definition | B: certified analytic decomposition | C: exact derived surface by definition |
| --- | --- | --- | --- |
| Charts | Finite spline parameter rectangles, interval evaluation and derivatives, regularity and injectivity proofs, and transitions across knots and periodic seams | Reuse analytic charts on each piece; prove their union covers the intended patch and bound each piece's deviation | Add charts for each derived law, with inverse/branch isolation and singular-domain handling; an evaluator alone is insufficient |
| Trims and boundaries | One shared 3D edge definition plus certified parameter traces on each incident face; complete trimmed-domain membership and orientation | Shared exact edges must join adjacent analytic pieces. More pieces mean more loops and real topology. Near-coincident independent boundaries are not a join | Retain exact intersections and branch selectors as trim definitions; certify finite order, closure and membership on the derived law |
| Predicates | Certified parametric membership and oriented-side queries. Use bounds as filters; retain exact or certified near-contact outcomes and distinguish the target from the candidate | Existing analytic predicates apply per piece; add complete cover and seam classification. They do not prove agreement with the source definition | Prove signs, equality, incidence and normals for the derived definition, including cross-kind analytic reductions and singular contacts |
| Intersections and booleans | Isolate all relevant roots and coincident regions in parameter domains, retain component certificates and parameter traces, and adapt all five boolean stages | Reuse supported analytic pair solvers, then reconcile meetings across every artificial edge. Face count increases arrangements and sewing work | Add derived/analytic and derived/derived solvers with complete branch coverage. Eliminating an offset to a polynomial can introduce unwanted branches and high degree |
| Measurement | Integrate the represented spline with certified derivatives. Enclose the ideal recipe properties by direct integration or the preceding ruling's tube and Jacobian proofs | Sum disjoint represented pieces with orientation; prove exact-target corrections and seam cancellation | Integrate the derived law over its exact trims. Exact geometry does not make quadrature free or remove property accuracy requirements |
| Meshing | Bound chord error from spline derivatives and share boundary samples; keep mesh-to-candidate error separate from candidate-to-target error | Reuse analytic meshing per piece, sharing all new boundaries; analytic creases remain actual geometry | Bound charts and derivatives of the derived law, cache only rounded display vertices, and cover every singular or trimmed region |
| Field measuring stick | Cover every trimmed spline span against the exact local definition with two-sided metric bounds, branch and reverse-coverage evidence. Existing analytic candidate coverage is insufficient | Measure every analytic piece and the union's coverage/gluing. Checking isolated pieces cannot detect a missing strip | Prove the represented law equals its exact target. The face's representation error is zero only with that proof; display and export still need their own deviation checks |
| STEP | Emit spline surfaces and shared curves with certified decimal rounding and trim compatibility. A degree-1 sweep may use `SURFACE_OF_LINEAR_EXTRUSION` | Emit analytic faces, all added topology and its correspondence; no hidden replacement of one smooth face with many creased faces | Use an eligible exact STEP form if available; otherwise fit only at export and certify the emitted solid. An exact in-memory law does not guarantee an eligible solid STEP representation |
| Names | Keep one face and its origins across knot refinement. Chart seams have no recipe edge identity; real trims/splits retain source origin sets | Every real split inherits the source origins. Artificial edges affect selection and adjacency and require an explicit topology/naming choice | Keep source origins through derived laws and analytic reductions. Branch/domain splits follow the same recipe split rule |
| Click cost | Local low-degree evaluation, reusable bounds and subdivision; fitting, complete intersections and proof may dominate. Exact extrusion reduces the first fit to one dimension | Cheap individual predicates, but potentially many faces, edge meetings, checks and export entities. Tightening the error can increase topology sharply | No fitting or model-error proof for an identity-preserving construction, but each derived law needs certified intersection, distance and integration paths; nested definitions can make these expensive |

Choose A because it provides a bounded, smooth candidate representation that can serve several feature families and has a corresponding STEP surface form. B changes the model's topology to avoid adding a surface kind, and its extra seams remain a cost after construction. C is mathematically sound and remains the exact fallback when all its consumers exist, but it postpones no intersection or trim obligation and does not deliver the first approximated face. No timing evidence here proves A faster than C. If A misses click targets, repair or reassess the measured bottleneck instead of weakening the proof.

An exact linear sweep of a certified spline curve is a useful specialization of A, not another feature definition. Preserve that degree-1 direction so the first implementation need not solve a two-dimensional fitting problem. It does not replace the general surface contract needed by fillets.

STEP supports knot and control-net data through [`b_spline_surface_with_knots`](https://downloads.steptools.com/docs/stp_aim/html/t_b_spline_surface_with_knots.html), and positive weights through [`rational_b_spline_surface`](https://downloads.steptools.com/docs/stp_aim/html/t_rational_b_spline_surface.html). STEP also defines [`offset_surface`](https://www.steptools.com/stds/stp_aim/html/t_offset_surface.html). Its existence must not be mistaken for permission to place it directly on an `ADVANCED_FACE`: [that entity's WR1](https://www.steptools.com/stds/stp_aim/html/t_advanced_face.html) restricts the support to elementary, swept or B-spline surfaces. Thus failure to fit the current six kinds and failure to export an exact definition are separate questions.

### Representation and topology contract for A

| Record | Required content and invariants |
| --- | --- |
| Exact definition | Immutable feature and input-definition references, exact `Real` dimensions and pose, operation law, branch/domain restrictions, orientation and origins. Record the whole-solid recipe target separately from the local face target. The dependency graph is acyclic |
| Spline surface | Positive integer tensor-product degrees `p,q`; exact rational nondecreasing knot vectors; rectangular control net of exact `Real` points; matching positive `Real` weights, with all weights one for a polynomial patch. Builders start with degrees 3 and 1 for the extrusion, and bicubic patches for general regular sheets. Prefer rational controls and weights when the constraints permit them. Higher degree consumes the same deterministic work allowance, never an unchecked special path |
| Domain | Finite parameter bounds, knot multiplicities and explicit closed/periodic transitions. For `n` control points in a direction, the expanded knot vector has length `n+p+1`, with active interval `[knot[p], knot[n]]`. Clamped endpoints repeat `p+1` times. Skip zero-width knot intervals; reject a zero-width active domain or inconsistent dimensions. Each rational denominator has a proved positive lower bound on the covered domain |
| Span cover | Ordered nonzero knot spans and, when needed, multiple joined charts under one logical surface. Store certified position, first- and second-derivative bounds, rank and branch information. De Casteljau/Bernstein bounds are proof tools, not samples. Prove injectivity on each regular chart and account for overlaps globally |
| Smooth joins | Shared exact boundary data, transition maps and matched first derivatives after reparameterization for internal joins of a smooth face. Positional and tangent continuity are proof obligations. Knots may subdivide computation without creating physical edges. A crease is a real boundary requiring the recipe's topology semantics |
| Trims | Shared 3D curves and vertices, oriented parameter traces on every incident chart, finite parameter ranges and certified correspondence. A general meeting curve may remain algebraic/implicit with root selectors; it need not be a spline. Do not independently fit a different boundary for each face |
| Local evidence | Candidate revision, exact face definition, covered trims, units, a proved worst-deviation enclosure `[L,U]`, feature-size bounds and proof dependencies. Record the conservative upper `U` as local tolerance. Neither a requested tolerance nor a sample maximum is this record |
| Body evidence | Valid solid, complete face/edge cover and gluing, outward material orientation, and a whole-result certificate against the recipe. Local records alone cannot create it |

These are mathematical data requirements, not a second certificate API. Coordinate their concrete types with the composable-error implementation and [the measuring-stick change](https://github.com/BohnBawerick/catallax/pull/30). This ruling does not change that change's scope or treat its current candidate evaluator as spline support.

A face has one represented geometry. Its predicates, intersections, mesh and STEP use that candidate consistently; a query against the ideal definition must say so. A later feature retains and constructs from the ideal feature definition, then certifies its candidate against the new ideal result. It must not offset last time's fitted control net and call that the exact target. A boolean on represented operands still owes the transfer proof in the preceding ruling. Consulting the original definition does not by itself prove that trimming the candidate selected the correct material.

For a trim shared with an analytic neighbour, prove that the common edge lies on both represented surfaces. A bound saying the edge is near both is insufficient. Fit with shared constraints, trim at the certified common intersection, or change the representation and reprove the affected faces. Preserve known construction incidence and tangency facts only when the candidate construction proves them. Normal closeness is not a new tangent fact. General fillets need a contact/continuity construction as well as a positional fit; if the chosen patch family cannot satisfy a required join, add a complete exact transition carrier or report unfinished coverage. Never hide a crack or a kink under the local tolerance.

Internal knots and chart cuts do not appear as selectable edges and do not split the origin set. A disconnected surviving domain or a physical boundary can split the face, and every resulting face inherits its origins under [RECIPE section 3.2](RECIPE.md#32-faces-carry-origin-sets). STEP may split the logical face into joined file faces where its representation requires it; retain an explicit mapping back to the logical face and certify those joins. Subdivision order and patch indices never become persistent recipe names.

### Consistency, caps and no growth

All three options obey PLAN's four rules. They differ in what must be proved at the representation boundary.

| Rule | A | B | C |
| --- | --- | --- | --- |
| Decide once | Key a bounded decision by immutable objects, domain and representation revision; all spline consumers read it | Every analytic piece touching the decision reads the same record | Exact construction facts answer exact contacts; any bounded downstream decision still gets the same shared record |
| Same becomes one object | Merge to one vertex or edge and shared trace evidence; revalidate all affected incidences | Share actual topology between pieces, including artificial seams | Share exact derived intersections and vertices, not repeated numerical evaluations |
| Record local tolerance | Store measured `U` against the exact face definition | Store each piece's `U` plus a complete union correspondence | Identity-proved exact faces have zero representation displacement; inherited bounded changes still need evidence |
| Check topology before publication | Check incidence proofs, trim coverage, shell orientation and intersections before the whole-result gate | Also check no missing, overlapping or duplicated tiles and no unrecorded crease | Check derived domains, singularities and selected components as well as the common shell invariants |

Merging is not transitive tolerance chaining. Every object absorbed into a representative needs evidence covering its displacement and the resulting incidences. A later exact query cannot contradict a recorded bounded choice about the same represented object. Rebuilding the candidate creates a new revision and requires proof that any reused decision still applies.

For kernel-made approximations in A or B, aim at `U <= 100 nm` and require `U <= 1 µm`. Also require `U <= l/10` and `U <= w/10`, using proved lower bounds on the shortest neighbouring edge length `l` and narrowest neighbouring face width `w`. Missing size evidence is unresolved, not permission to omit the cap. Computational knot spans are not neighbouring physical features; B's real new edges and faces do count. If a cap fails, refine, split/refit or substitute exact geometry; never clamp `U`. C introduces no local representation displacement when equality is proved, so those caps hold at zero. Its fitted exports remain a separate representation.

Imports have the PLAN limit of 10 µm unless the model records an explicit file-specific acceptance of more. Any chosen carrier preserves that acceptance and distinguishes the imported geometric definition from stated source uncertainty. The import allowance does not grant a kernel-made approximation 10 µm. List each imported face and its bound in the complete result certificate; hand-offs retain that scope and add their own proved output error.

The whole-result certificate checks the local policy over a complete valid solid. It does not conjoin the face caps with the superseded universal sub-1 nm bound:

```text
exact region: equality proof against the retained definition
approximated/imported face: proved U_i within its recorded local cap and feature-size caps
bounded decision: substituted gap < 1 nm, with its affected scope and consequences proved
whole body: complete exact/exception cover, shared incidence and material orientation
output region: declared accuracy and composed model-to-definition plus output error
```

A certified 50 nm fitted face can meet the 100 nm target and pass the body certificate. Its cap record alone is insufficient: coverage, gluing, correspondence and orientation must also pass. An unlisted changed face, a 2 nm sameness decision or a reversed shell still fails. Exact faces do not inherit a neighbour's allowance. A uniform signed-distance bound may be derived when useful, but it is not another sub-nanometre gate. The caps and the 1 nm decision tolerance serve the two permitted approximation cases in PLAN section 3.2; they are not tolerance-by-purpose tiers for geometry creation.

No growth means remeasure the current candidate against the current exact definition. For A, refit from that definition when needed; for B, rebuild the analytic cover; for C, compose exact laws with their branch and domain proofs. Restricting a certified face to a smaller valid domain can reuse an upper bound, but changed trims and newly adjacent short features can require a tighter one. A radius-addition identity may simplify an offset chain only on a proved regular branch; medial events and changing nearest witnesses cannot be ignored. Resetting a local record never discards a whole-result dependency. Only a new direct proof against the recipe, or the approved transfer proof, can replace the inherited bound.

Failed proof budgets and unsupported carriers remain unfinished coverage. A geometrically self-touching result that violates the recipe's solid rules remains a correct impossible-geometry refusal. An inward offset proved to remove everything returns the exact empty result. A singular chart alone proves neither case.

### First feature and build order

Land the first tolerance-core piece with a real approximated offset face, its local field proof and the whole-result gate. Foundation-only commits may support review, but do not count as the feature's completion. The stages below describe one coordinated implementation across consumers, not permission to publish a partly adapted surface kind.

1. Use the approved carrier and the amended local-certificate policy. Retain the existing exact six-kind paths and inspect the final measuring-stick interface when it lands. Define the spline/trim/certificate schema once across shapes, topology, intersections, measurement, field, mesh, features and export.
2. Use a closed right elliptic cylinder with `a=2 mm`, `b=1 mm`, `0 <= z <= 4 mm`, and inward offset `d=-1/4 mm`. Build the source through the public recipe path: extrude a radius-2 circle in a plane tilted 60 degrees about the x axis, centred at `z=-2`, by 8 mm along z; intersect it with the box `[-3,3] x [-3,3] x [0,4]`. Its horizontal cross-section has semi-axes 2 and 1. This needs no new ellipse-profile syntax, but any missing source-path coverage must land with the feature. The offset's planar caps move to `z=1/4` and `z=15/4`. The ellipse's minimum curvature radius is `b²/a=1/2 mm`, so its inward cross-section stays regular and convex. Erosion of this Cartesian product is the eroded cross-section times the shortened interval. It needs the new side carrier and planar caps, with no outward edge-round construction hidden offstage. The earlier outward `d=1` example remains the regression against adding ellipse radii; it is not a claim that its entire capped-solid offset needs no round faces.
3. Construct the exact inward face definition and a closed cubic cross-section fit with exact periodic position and tangent joins, then extrude it with degree 1. Prove convexity, injectivity and finite-domain coverage. Share its end curves exactly with the planar caps. Keep the source side/cap origins. Implement chart inversion, membership, oriented queries, supported intersections and finite trim traces; replace the assumption that every surface has one usable global `RealPoly`. Preserve analytic fast paths. Unhandled pairs must give specific coverage refusals, never an invented analytic type.
4. Extend the field's candidate-cover consumer to this spline representation, using the separately owned measuring stick. Certify two-sided distance to the exact local definition over every retained span and boundary, including the nearest-witness branch. A distance-to-source residual needs the regular normal-neighbourhood/correspondence proof before it becomes offset-face distance. Include the planar caps' changed trim domains: their exact plane supports do not prove that their fitted perimeters are exact. Require the local cap and feature-size checks. A `converged=false` or a coarse supporting-domain enclosure is no success unless its actual upper bound meets every applicable gate; samples only supplement the proof.
5. Close the whole-solid proof for this fixture. For a compact convex body with interior, `s(x) = sup over unit u of (x dot u - h(u))`. Outside, the nearest point supplies a supporting normal attaining the distance; inside, the smallest supporting-plane gap is the radius of the largest ball about `x` contained in the body. A uniform support-function bound of `epsilon` therefore bounds the oriented signed-distance difference by `epsilon`. Prove the fitted and ideal cross-sections convex and bound their support difference over every direction; the shared exact axial interval preserves that bound for the product solids. This can establish a useful whole-solid distance summary alongside the required local records and material proof. Fit to the 100 nm target within the 1 µm hard cap and feature-size limits; do not force `epsilon < 1 nm`. Choose export accuracy separately and report its composed bound. Stop model refinement when the local policy and complete body proof pass, unless a tighter requested output or property needs it.
6. Complete mesh, properties, STEP and application exchange for the new face. Use the spline derivative bounds for mesh coverage and certified target-property integration or corrections. Emit spline/sweep surfaces with common end curves, then certify the rounded file. Add a separately reviewed versioned kind/tag and seam/trim contract where required; the legacy exchange forbids pcurves and seam curves and has fixed surface tags. Do not masquerade as an ellipse or silently widen that exchange. A returned body must be usable by all these consumers, with property work deferred until requested as PLAN requires.
7. Exercise a radius edit, a certified planar trim/boolean and a second offset from the retained original definitions. Recheck origins, shared boundaries, local caps and operation transfer. Include decisions just below/at/above 1 nm, an accepted 50 nm fit, fits at/above the inclusive 1 µm cap, the import limit and recorded override, a short neighbouring feature, wrong normal/branch, missing reverse coverage, a duplicated or missing span, a forged seam fact and a stale certificate. Test true empty erosion, invalid self-touching results and insufficient-work refusals separately. Read STEP independently, including decimal surface/edge compatibility. Repeat deterministic evidence at one and two workers locally and four in CI; compare cold/warm caches without changing logical work allowances.

The first feature is complete only when the application/verification path returns this valid solid with a genuine spline side face, local field evidence, a complete whole-result certificate listing its bounded faces, display, stable origins after an edit, requested properties and certified STEP. Inspecting or measuring a pre-existing analytic face against an offset definition does not satisfy that gate.

Next, add outward elliptic-cylinder offsets with their edge pipes and corner strata, then general rolling-ball fillets from the selected offset-contact path and envelope in RECIPE section 4.5. Keep exact cylinders, spheres and tori wherever the construction reduces to them. Fit only the remaining regular sheets, preserving shared contact curves, continuity constraints, ends, corners and consumed-face history. Reuse the spline consumers and proof records; extend the exact definition and field metric for each new family. General surface pairs, singular envelopes and contact changes still need new proofs. A successful ellipse fit proves none of them automatically.

Geometric blends follow their specified exact laws. Field smooth-min blends remain the separate field-only feature in PLAN section 3.6; do not rename one as a rolling-ball fillet. Their conversion, when requested, uses the same spline carrier only after a complete definition, trimming, topology and deviation proof exists. Helices or transport laws outside algebraic `Real` still require the previously identified certified-function extension.

Measure end-to-end rebuild, first complete display, local certification, whole-result certification, requested measurement and STEP separately as well as together. Retain the architecture targets of 1,000 ms for feature builds, 500 ms for booleans and 2,000 ms for background comparison including requested output work. Cache immutable definitions, knot-span bounds and verified certificates by revision and domain; rebuild only affected spans. Initial fitting and certification stay on the publication path. A background preview or incomplete proof cannot be reported as a successful solid. Timings in a shared environment are diagnostic; use warmed, pinned, matched runs for a performance claim.

### Approved PRD change

The owner approved the replacement of [PRD section 3.1 item 1](PRD.md#31-agreement), and this change applies it. The former sentences "Both kernels hold these exactly by their definition. The approximation appears only in the exported file." no longer prohibit a fitted B-rep carrier. Both kernels still retain exact definitions. Each approximated/imported face carries its measured local tolerance within the applicable cap, and the body certificate covers the exact regions and all exceptions.

The section's introductory agreement sentence now refers to those scoped bounds rather than a universal 1 nm difference. Item 2 retains the strict 1 nm decision-gap rule. Measurement, comparison and STEP use the amended composable-error policy above; an existing exchange promising a stronger bound remains unchanged until a versioned contract explicitly adopts the new policy. [PLAN section 14 item 3](PLAN.md#14-proposed-prd-changes) records this decision.

### Observed baseline and evidence for this ruling

Source observations below are at `046f426b8f42a08a32a7ae88f579dd32482b273e`. Tests were read for their assertions, not executed for this docs-only ruling. Design requirements above must not be read as implemented behaviour.

| Observed behaviour | Inspected code and tests |
| --- | --- |
| `Surface` has six analytic kinds, analytic `implicit`, `gradient`, `domain`, `evaluate` and `derivatives` methods. No spline or derived offset carrier exists | [surface.rs](../crates/catallax-shapes/src/surface.rs); [algebraic.rs](../crates/catallax-shapes/tests/algebraic.rs), `algebraic_dimensions_preserve_incidence_gradient_and_implicit_scale`, checks exact equations and gradients |
| A face stores one `SurfaceId`, sense and loops. The topology checker reads incidence facts structurally; it does not reprove geometric incidence | [entity.rs](../crates/catallax-topology/src/entity.rs); [check.rs](../crates/catallax-topology/src/check.rs); [geometry.rs](../crates/catallax-topology/tests/geometry.rs), `intersection_edges_and_meeting_vertices_need_no_coordinates` and `check_uses_factbook_closure_not_only_direct_records` |
| An intersection curve names two surfaces and a regular witness. Finite curve order and branch handling are separate consumers | [curve.rs](../crates/catallax-shapes/src/curve.rs), `IntersectionCurve`; [chart.rs](../crates/catallax-intersect/src/chart.rs); [conic_sections.rs](../crates/catallax-intersect/tests/conic_sections.rs), `open_conics_have_ordered_finite_traces_and_correct_contacts`; [boolean geometry.rs](../crates/catallax-boolean/src/geometry.rs), `on_arc`, refuses unresolved ordering |
| Measurement has a dedicated elliptic-cylinder integrand, and the cross-crate fixture checks its mass properties. This is not general spline integration | [measurement surface.rs](../crates/catallax-measure/src/surface.rs), `compile`; [STEP algebraic.rs](../crates/catallax-step/tests/algebraic.rs), `oblique_circle_extrusion_has_exact_measurements_and_step_geometry` |
| Display `Chart::new` explicitly refuses elliptic cylinders at this baseline. Existing seam/pole tests cover block, circular cylinder, sphere, cone/frustum and torus | [mesh chart.rs](../crates/catallax-mesh/src/chart.rs); [mesh tests.rs](../crates/catallax-mesh/src/tests.rs), `primitives_close_at_seams_poles_and_apices_with_bounded_chords`. Adding a shapes evaluator has not by itself completed display coverage |
| STEP emits an elliptic cylinder as `SURFACE_OF_LINEAR_EXTRUSION`; its fitter produces degree-seven curve spans, not surface patches | [export.rs](../crates/catallax-step/src/export.rs), `surface`; [fit.rs](../crates/catallax-step/src/fit.rs), `fit`; [algebraic.rs](../crates/catallax-step/tests/algebraic.rs), the same extrusion test reads the emitted entity and ellipse data |
| Field query types distinguish true distance, residuals and geometric distance bounds. `face_query_checked` checks named-face support and surviving boundary, not an arbitrary spline-face supremum | [field lib.rs](../crates/catallax-field/src/lib.rs), `Value`, `DistanceBound`, `face_query_checked`; [metric.rs](../crates/catallax-field/src/metric.rs). The [measuring-stick proposal](https://github.com/BohnBawerick/catallax/pull/30), read on 2026-10-23, reports analytic-chart whole-face covers and explicitly excludes fillet/offset construction |
| Recipe evaluation falls through to `features.unsupported-operation` for unimplemented operations. Origins are ordered by canonical recipe content | [evaluate.rs](../crates/catallax-features/src/evaluate.rs); [names.rs](../crates/catallax-recipe/src/names.rs), `Origin::key`. No offset-construction success is inferred from the presence of recipe definitions |

The earlier [error-policy decision](https://github.com/BohnBawerick/catallax/pull/14) supplied the composition and publication obligations. The owner-approved amendment above replaces its universal sub-1 nm acceptance gate with complete exact/exception coverage and local caps. The measuring-stick work supplies local evidence; it does not by itself establish that complete body certificate.

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
question/manifest exchange. Display chord bounds do not establish the model's
geometric certificate.

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

The table describes observed inherited behaviour, not completion of the new PRD. The amended composable geometric error decision above defines conversion, composition, the strict 1 nm decision limit and the local approximation caps. Enforcement remains implementation work; this publication changes no numerical code or outcome API.

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

Both evaluators use the same operation definitions and prove exact agreement on exact regions. PRD section 3.1 permits bounded decisions below 1 nm and approximated/imported faces within their recorded local tolerances and caps. The B-rep retains the exact definition beside a fitted carrier; STEP carries separate output error. The complete certificate lists every exception.

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
