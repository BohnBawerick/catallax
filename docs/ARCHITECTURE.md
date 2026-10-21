# Catallax architecture

Catallax's master model is one recipe evaluated by an exact B-rep kernel and a signed distance field kernel. Each operation shares one mathematical definition. [PRD.md](PRD.md) defines that direction; [DESIGN.md](DESIGN.md) governs the B-rep implementation.

The crates below are the implemented B-rep half, exported without source history into this public repository. The recipe language, field kernel, cross-kernel checks and conversions are not implemented. Their boundaries must be defined before adding crates for them. [PLAN.md](PLAN.md) orders that work.

The target is reasonably good, comparable to OpenCASCADE, SMLib and the Rust kernels. The former 1-in-10,000 programme and grid-cell product gates are retired. Mathematical regression tests remain; speed per click matters for the interactive CAD tool.

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
Future recipe evaluation calls these existing B-rep operations and a separate field evaluator. The two kernels do not convert results after each operation. Conversions happen only when a feature or output needs them.

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
crates/catallax-cli/             front end, binaries `catallax` and `exact-brep`
docs/                      this file, DESIGN.md, JUDGE-INTERFACE.md, INTERFACE-CHANGES.md, pieces/
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

Treating `NotExact` as coincident is the single rule that turns an unfinished proof into a bounded, visible tolerance. No caller may instead pick a side from the cached approximation.

The table describes inherited behaviour, not completion of the new PRD. The required bounded outcome is limited to a proven geometric gap smaller than 1 nm, or 0.000001 mm. Scalar and angle marks need a proven distance interpretation and composition across the result. Enforcing that cap is still open; this publication changes no numerical code or outcome API.

Only impossible geometry is an acceptable final refusal. Unsupported pairs, ambiguous classifications, exhausted budgets, deadline refusals and missed measurement accuracy are unfinished implementation work. Preserve their honest stable codes until the underlying coverage is repaired.

### Operation level: `Answered<T>`

Every operation (an intersection, an imprint, a boolean, a measurement, an export) returns `catallax_numbers::Answered<T>`, which is `Result<Answer<T>, Refusal>`.
`Answer<T>` holds the value and a `MarkLog`: every `NotExact` met on the way, each turned into a `Mark` with its reason code, its bound and the kind of bound (a distance, an angle, or the raw quantity).
An answer with an empty mark log is exact. An answer with marks is "marked not exact with proven bounds", and the marks travel with the result: `Body` stores the marks of the operation that produced it, merged with the marks of its operands.
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

Open item: the CLI's 50-second `numbers.deadline` cutoff can answer or refuse the same question depending on machine load. That conflicts with the same-answer rule across machines and thread counts. Resolve the policy and exchange contract separately; no deadline code changes here.

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

Future persistent face and edge names identify the creating operation and its inputs, never the traversal order or an arena index. Store-local `FaceId` and `EdgeId` handles remain useful during evaluation but are not rebuild-stable recipe references. Define splitting, merging and feature-edit naming rules with the recipe language.

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
| 4 | 5, 6, 7 | `SurfacePair`, `Component`, `Contact`, `CurveSurface`, `CurveHit`, `Cyclic`, `CurveForm`, `CurveTrace`, `TracePiece`, `intersect_surfaces`, `intersect_curve_surface`, `curve_form`, `on_curve`, `same_curve`, `order_on_curve`, `cyclic_order`, `tangent_at`, `arc_point`, `trace`, `supports` |
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
