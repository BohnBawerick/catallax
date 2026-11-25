Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Exact and tolerant in one kernel: how Catallax works today and the coexistence options

Repo: Catallax repository at bb2f8375. Read-only.
Labels: [code] = read in source; [doc] = read in repo docs; [notes] = earlier scratchpad notes or grill pack;
[memory] = literature from memory, not re-read this session; **Inference** = my reasoning.

## 1. Today: how exact and numerical bodies interact

### 1.1 Two types, two crates' worth of contracts

- Exact body: `catallax_topology::Model` + `BodyId`. Edges hold only the 3D curve (`crates/catallax-topology/src/entity.rs:44-62`) [notes]. Checked by `Model::check` (`crates/catallax-app/src/lib.rs:482-489`) [code].
- Numerical body: `catallax_topology::numerical::NumericalBody` (`crates/catallax-topology/src/numerical.rs:117-140`) [code]. It is a separate definition (vertices, edges, trims, native trims) with ONE tolerance for the whole body: `NumericalOptions { tolerance_mm }`, default 0.01 mm, "not a promised global deviation bound" (`numerical.rs:100-115`) [code].
- The doc comment states there is "no conversion to an exact evaluation checkpoint", and a `compile_fail` test forbids `NumericalBody -> Answer<NumericalBody>` (`numerical.rs:117-140`). Same fence on `NumericalEvidence -> BoundedEvidence` (`crates/catallax-numbers/src/numerical.rs:52-62`) [code].
- So there are three result tiers: exact `Answer` with empty marks; bounded `Answer` with a mark log (`crates/catallax-numbers/src/outcome.rs:205-217`); numerical `NumericalEvidence` (named sample checks, no bound between samples) [code].

### 1.2 The one-way bridge: exact -> numerical

- `catallax_features::numerical_embedding::definition_from_exact_model` (`crates/catallax-features/src/numerical_embedding.rs:907-940`) [code]: takes an exact `Model` body, keeps native carriers and curves, fits traced intersection edges as cubics, fits u/v curves, sets vertices to 128-bit midpoints, converts the tolerance with `to_f64_nearest` (`:927`), then `build_with_native_trims` must admit the result.
- Callers [code, grep]: numerical fillet (`evaluate_numerical_fillet.rs:127, 228`), rolling-ball bodies (`fillet_rolling_body.rs:443`, `fillet_rolling_vertex_body.rs:334`, `fillet_native.rs:1970`), helical/general sweeps (`crates/catallax-app/src/numerical_sweeps.rs:120`), numerical STEP and loft STEP in the CLI (`numerical_step.rs:210, 309`, `numerical_loft_step.rs:108`).
- No numerical -> exact path exists (by design, see 1.1).

### 1.3 What happens when an operation meets one of each

| Case | Today | Evidence |
|---|---|---|
| Any recipe step after a numerical fillet that reads that body (Boolean, extrude onto its face, second fillet, pattern) | **Refused**: `features.numerical-continuation`, "its native numerical continuation has no admitted construction for this operation" | `crates/catallax-features/src/evaluate_numerical_fillet.rs:47-66`, called at `evaluate.rs:890-892` [code] |
| Recipe outcome | Either `EvaluationOutcome::Exact(checkpoint)` or `Numerical(body)`, never mixed | `evaluate_numerical_fillet.rs:7-14`; app `lib.rs:455-480` [code] |
| Boolean on a fitted sweep/loft | Fitted rational charts are held inside an exact `Model`; the **exact intersector** runs on them ("Boolean operations still use the fitted definitions and retain the existing intersector's coverage refusals"); the whole result is embedded to numerical at the end | `numerical_sweeps.rs:15-20, 108-140`; `docs/RECIPE.md:992` [code, doc] |
| Numerical Boolean module | "Numerical material selection after an exact finite intersection and imprint"; output is a numerical body via callback, "never an exact Boolean answer" | `crates/catallax-boolean/src/numerical/loft.rs:1-3` [code] |
| Fillet on an imported part | Recipe import goes through strict `catallax_import::apply_source` into the exact `Model` (`evaluate.rs:916-924`); then the fillet tries exact, refunds, embeds, builds numerically. The measured case exhausted 98.4M in import and never reached embedding | `crates/catallax-import/src/repair.rs:79-85`; `docs/phase5-work-cases.md:86-95` [code, doc] |
| Tolerant import in the app opener | Per body `Tier` (ValidSolid, RepairedSolid, TolerantSolid, SurfaceBody) and an `approximate` flag; storage `Certified{model}`, `Numerical(StepImport)` or `Tolerant` | `crates/catallax-import/src/tolerant/opened.rs:19-53` [code] |

Summary: there is **no promotion and no mixed body**. Exact is converted wholesale to numerical when a numerical step is needed (one-way bridge), and anything after that is refused, except sweeps and lofts, which cheat by keeping fitted charts as exact rationals inside the exact Model and running exact arithmetic on them.

### 1.4 f64 <-> rational

- f64 -> exact: `Rational::from_f64_exact` (`crates/catallax-numbers/src/rational.rs:233-257`), lossless (power-of-two denominator) [code]. Used to store every fitted control point, knot and weight (`numerical_embedding.rs:57-58`, `fillet_numerical.rs:51`) [notes].
- exact -> f64: `Rational::to_f64_nearest` (`rational.rs:259-262`), ties to even [code].
- Pattern: "Binary64 arithmetic only proposes fits. Numerical admission must recheck every edge, trim and loop at 128 bits" (`numerical_embedding.rs:5-8`) [notes]. So even numerical bodies store geometry as exact rationals; only the proposal and the 128-bit ball checks are inexact.
- **Inference:** this is why exact operations on fitted geometry are slow: they run rational arithmetic on 53-bit-mantissa dyadic numbers of high degree, so bit length explodes.

### 1.5 Label propagation

- Per body today: `NumericalEvidence` per `NumericalBody` with one tolerance; `BodyOutcome::approximate()` per imported body (`opened.rs:46-53`) [code].
- On paper, per face: the owner's 23 Oct amendment: "Certify the whole result as exact except for its listed approximated or imported faces and recorded bounded decisions. Each exception carries a proved local tolerance" (`docs/ARCHITECTURE.md:291`). Exact faces have zero displacement (`:541`); a planned test is "an exact region beside a valid 50 nm approximated face without granting the exact region 50 nm" (`:459`); "merging is not transitive tolerance chaining" (`:544`) [doc]. Enforcement is open (`docs/DESIGN.md:28`) [doc].
- Grill pack D7 (per face vs per part) is still open [notes].
- Caps: 1 nm modelling sameness (`crates/catallax-numbers/src/certificate.rs:20-38`), local caps 1e-4 / 1e-3 / 1e-2 mm (`certificate.rs:262-277`) [code].

## 2. Cost profile

Exact is cheap:
- Box, all 12 edges filleted exact: 9,467 to 13,195 work (`docs/phase5-work-cases.md:489`).
- Rolling-ball rim strips that land on native cylinders and tori: ~6M, flat in radius (`phase5-work-cases.md:463-474`).
- Certified-only imports: median 83,078 open work, 3% export share (grill pack, census CSV) [notes].
- Plane/quadric closed forms: two planes or plane+perpendicular cylinder give exact cylinder/torus blends (`docs/RECIPE.md:589`).

Exact is slow or fails:
- One fillet-end vertex, exact: 13.6M (`ROADMAP.md:32`; `docs/checkpoints/checkpoint-1.md:272`).
- Exact STEP export of curved faces: 8.5 s at only 807,063 work, clock in exact rational Taylor controls (`checkpoint-1.md:274`). Work is not time: 95 to 448,000 units/ms [notes].
- Exact intersector on fitted geometry: M8 later cuts on a loft body, 80M to 100M each, most `numbers.operation-budget`; others `intersect.unsupported-pair`, `boolean.parametric-query` (`phase5-work-cases.md:125-139`). Real housing stops at `cut` (`:116-120`).
- Point containment on a four-contact fan: 27M to 53M per query (`phase5-work-cases.md:177-188`).
- Numerical admission is the cost in numerical builds: loft fillet 96.7M admission vs 0.63M geometry; P5 embedding+admission 98% (`phase5-work-cases.md:158, 308-312`) [notes].
- Numerical STEP re-admission: rim 21.9M to 22.6M, 41.8M near limit (`phase5-work-cases.md:484-488`).
- Numerical imports: median 1,663,900 open work (20x certified) [notes].
- Near-tangent: tangent contacts refuse by name (ROADMAP intersection row, `ROADMAP.md:77`, `:295-297`).

Pattern (**inference**): exact pays off when inputs are low-degree with small rationals (planes, quadrics with rational parameters, closed-form blends). It loses when inputs are fitted splines carried as dyadic rationals, near tangency, or degree-raising constructions (STEP Taylor controls).

## 3. Literature

- Jackson 1995 (Parasolid tolerant modelling) [read, scratchpad/jackson/j.txt]: "allows ... to mix geometries of different tolerances in the same model" (l.55-60). "Tolerances are associated with topology, not geometry. Points, curves and surfaces are treated as 'exact'" (l.336-340). Edge at default tolerance = one 3D curve; tolerant edge = pcurves in a tube (l.176-191). Vertex tol >= edge tol >= face tol (l.196-198). Coincidence uses the sum of tolerances (l.340-345). Jackson rejects exact rational arithmetic for "build-up of complexity of the numerical values" (l.46-49).
- Parasolid docs: exact entities behave as half session precision, 5e-6 mm; "enquiries on exact geometry are usually faster than enquiries on approximate SP-curves"; `PK_EDGE_optimise` tightens tolerant edges [notes: tolerance-numbers-research.md:41, t-smarter-research.md:19-60].
- ACIS: EDGE (3D primary) and TEDGE (pcurve primary, lazy 3D) coexist; tweak re-intersects tolerant edges to exact [notes].
- CGAL / Shewchuk: filtered predicates (float, then interval, then exact) over float inputs [memory]. Gives consistent decisions about the stored data, not about design intent.
- Hu, Patrikalakis, Maekawa 1996, interval solid modelling: geometry as intervals; exact = zero width [memory].
- Fortune 1997 / Milenkovic: exact compute then round with repair, polyhedra only [memory].
- Catallax already has the filter pattern: binary64 intervals first, exact fallback, in mesh and field (`docs/CAPABILITIES.md:159-162`) [notes].

## 4. Possibilities

### (a) Exact as tolerance zero in one topology type (Parasolid form)
How: one `Model`. Vertex/edge/face gain `tolerance: Option<Rational>`; `None` = exact (acts as a floor, say 5e-6 mm like Parasolid, or 1 nm). Exact builders write `None`; fitted, traced and imported entities write measured gap x1.1. Geometry stays stored exact (it already is, even for fits). Coincidence: if both exact, today's exact predicate; else distance vs t_a + t_b.
Pros: mixed operations become possible (Boolean exact with imported, fillet on import, extrude onto fitted face). Label per entity for free; matches the owner's 23 Oct per-face exception ruling. Kills `features.numerical-continuation`. One checker.
Cons: every exact algorithm needs a tolerant branch for tolerant inputs; tolerance spreads to vertices and edges of exact faces that touch a tolerant edge (Jackson's vertex >= edge rule); if tolerant geometry still goes through exact rational arithmetic, M8-style blowups remain.
Must change: `entity.rs` fields; `Model::check` absorbs `numerical.rs` checks; `NumericalBody` becomes a view or is removed; Boolean classification and imprint take tolerances; STEP writer mixes exact and tolerant entities in one file; measure tracks per-face U.

### (b) Two worlds, one-way bridge (today, extended)
How: exact until first contact with fitted or imported geometry, then `definition_from_exact_model` converts the whole body; every later op runs in the tolerant world.
Pros: exact code untouched; already built; clear boundary.
Cons: one fitted face makes the whole part tolerant forever (label overstates, drawings mark exact holes approximate); every op must exist twice (exact Boolean and tolerant Boolean, exact and numerical STEP writers: already true); the tolerant world today implements almost nothing after a fillet; results depend on step order; edits flip a body between worlds.
Must change: a full tolerant operation set (Boolean, fillet, extrude, pattern on `NumericalBody`); tolerant checker; per-body label.

### (c) Exact predicates inside a tolerant kernel (CGAL/Shewchuk style)
How: geometry f64 (or f64 copies of stored rationals) with per-entity tolerance; constructions in f64; each decision (sign, ordering, inside/outside, "distance <= t_a + t_b") through a filter: f64, then interval, then exact on the f64 inputs.
Pros: fast constructions; decisions consistent, so no contradictory topology from one comparison; same answer on every machine; reuses `interval.rs` and the mesh filter pattern.
Cons: exactness is about stored data, not intent: a designed coincidence (two coaxial cylinders built separately) is only found by tolerance, not by proof; closed-form exact blends degrade unless kept separately; near-tolerance comparisons fall to exact fallback often (cost in the hard middle).
Must change: intersectors and builders gain f64 paths; decisions route through filtered predicates; exact rational storage kept or dropped per entity.

### (d) Drop exact except where cheaper (fixed rule or predictor)
How: rule by geometry kind: planes, lines, circles, rational-parameter quadrics and known closed forms run exact; splines, fitted, near-tangent (normal angle test) go tolerant directly, without the exact attempt and refund. Optional per-step memory as predictor (grill Q1 = C).
Pros: removes the 13.6M vertex, the M8 100M cuts and 8.5 s STEP class; keeps the 9k box fillet; no doomed exact attempts.
Cons: needs (a) anyway to put exact and tolerant entities in one body; the rule boundary is where mixed-mode bugs sit; an edit that crosses the boundary (radius, near-tangent) can change the result kind; a predictor must be calibrated in time, not work.
Must change: path choice table in features; tangency flag; drop exact-first in import for spline-heavy files.

### (e) Others
- e1 Exact definition plus tolerant representation (ARCHITECTURE 305-324): keep the recipe solid R as the reference, the body B tolerant; tolerance = proved distance B to R; exact faces U = 0. Fits tolerance-first if the "proved" part is relaxed to "measured x1.1". Gives a path to later certify.
- e2 Tighten down (Parasolid `PK_EDGE_optimise`, ACIS tweak): after a tolerant op, try to re-intersect tolerant edges between exact surfaces exactly and drop the tolerance. Keeps exact islands from shrinking over a long recipe.
- e3 Lazy upgrade: tolerant first for speed, exact proof later on request (export, measure). `VISION.md:42` forbids durable proof jobs, so only inside one request [notes].
- e4 Interval solid modelling: tolerance as interval width. Uniform but heavy; no shipping kernel uses it [memory].

## 5. Recommendation (inference)

Adopt (a) in the Parasolid form, with (d) as the fixed rule for which constructions produce exact entities and (c) for how decisions on tolerant entities are made. Use e2 to keep exact regions from eroding. Treat (b) as the migration path, not the end state.

Reasons:
1. Jackson's split matches Catallax storage already: geometry is stored exact even for fits (`from_f64_exact`); only topology needs a tolerance field. So (a) is a topology change, not a geometry rewrite.
2. The owner already ruled per-face exceptions in an exact result (ARCHITECTURE.md:291). (a) is that ruling with measured rather than proved tolerances.
3. (b) costs the most long-term: double operations forever, and the per-body label overstates. Today's `numerical-continuation` refusal is (b)'s main failure, and it is live.
4. Speed: the slow cases are exact arithmetic on fitted or near-tangent data. A tolerant entity must be computed in f64 with filtered decisions, never fed to the exact intersector. Exact stays only for exact-exact pairs, which is where it is fast.

Rules to state up front:
- An exact entity never takes a tolerance except a vertex or edge forced by a tolerant neighbour (Jackson vertex >= edge >= face). Faces keep exact surfaces.
- Exact-exact decisions stay exact; any decision with a tolerant participant is "distance vs sum of tolerances" through a filtered predicate.
- Label is per entity; body and property labels are derived (a mass is approximate only if a tolerant face contributes).
- Same recipe, same machine, same path: the rule is by geometry kind, not timing, so results do not depend on load.

Risks: mixed-mode bugs at exact/tolerant borders (an exact vertex at the end of a tolerant edge); tolerance growth through chained Booleans; two STEP writer paths in one file; tests must cover the "exact region beside a tolerant face does not inherit the tolerance" case (ARCHITECTURE.md:459).

Open: whether exact entities use tolerance None (true exact) or a floor (Parasolid 5e-6 mm). None keeps the existing proofs; a floor simplifies the tolerant branch. Measure the M8 cuts under both before choosing.
