Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Choosing exact or fitted per operation: predictor research

Repo: Catallax repository at bb2f8375. Read-only.
Labels: [read] = read this session in source or repo docs; [snippet] = seen only as a search-result
summary; [unread/memory] = from memory, not checked this session; **Inference** = my reasoning.
Outside consult: `grok -p` (skeptic) was asked blind and hit its usage limit. No outside view.

## 0. What Catallax already does (the starting point)

- Exact-first with refund already exists for fillets. `evaluate.rs:1490-1515` runs the exact builders
  inside `catallax_numbers::refund_on_refusal`; on refusal the rolling-ball construction runs with the
  allowance the attempt found [read].
- `refund_on_refusal` (`crates/catallax-numbers/src/operation.rs:192-233`) restores `used` and clears an
  exhaustion the attempt latched. Its doc says: "The arithmetic of a refused attempt is real machine
  time, so the caller bounds how often it makes one." It also marks the run-ahead tape unrecordable [read].
- **Inference:** so the attempt is uncapped. A doomed exact attempt can burn up to the whole remaining
  100M of machine time, then disappear from the reported work count. Today's "exact-first" is design (c)
  with k = infinity and no record of the waste.
- The work meter charges decider arithmetic by limb count (ball add 4L, multiply 8L^2;
  `decider.rs:1659`, `limbs()` at `:1834`) [read]. But `crates/catallax-step/src/fit.rs`
  (`taylor_controls`) has no `charge` call at all [read]. That is why exact STEP export of curved faces
  took 8.5 s on 807,063 work (`docs/checkpoints/checkpoint-1.md:275`) [read]. Prior notes put the
  observed spread at 95 to 448,000 work units per ms [notes, not re-verified].
- **Inference:** any predictor or cap counted in work units is blind wherever the meter is blind. Meter
  fidelity is a precondition, not a detail.

## 1. Literature and practice

### Algorithm selection and hardness models
- Rice 1976, "The algorithm selection problem", Adv. Comput. 15: problem space, feature space, algorithm
  space, performance mapping, selection mapping S(x) [snippet, via Kotthoff].
- Kotthoff 2014, AI Magazine 35(3) 48-60, survey [snippet]. Portfolios correspond to Rice's algorithm space.
- Leyton-Brown, Nudelman, Shoham 2009 JACM (EHM methodology); Hutter, Xu, Hoos, Leyton-Brown 2014 AIJ 206
  "Algorithm runtime prediction" [snippet]. Forward selection found 5 to 10 features often predict as well
  as all; tree/random-forest models did best. They train on thousands of instances per domain [snippet
  for features, unread/memory for instance counts].
- **Inference:** EHMs work on large instance sets with log-runtime targets and accept large per-instance
  errors. Catallax has about 40 measured rows, dominated by two parts (owner rim, loft). A learned model
  would fit the test parts. That is exactly the risk the owner named.

### Filters, degeneracy, lazy exact
- CGAL filtered predicates and Lazy_exact_nt: interval filter first, exact only on filter failure; the
  lazy kernel keeps a DAG and replays it exactly when needed (Fabri and Pion, arXiv cs/0608063) [snippet].
- Pion 2004 slides "Robustness issues in CGAL": counts of static and dynamic filter failures on Delaunay of
  1e5 points on a 30-bit grid; "filters are fragile: try to avoid bad cases in algorithms" [snippet,
  numbers not seen].
- Failure rate depends on the algorithm's predicate choice, not just input (Voronoi via envelopes,
  arXiv 0906.2760) [snippet].
- Shewchuk 1997 adaptive precision predicates: stages of growing precision, cost paid only as far as the
  sign needs [unread/memory].
- Attene 2020 / Cherchi et al. 2020 and 2022 indirect predicates: implicit points (LPI, TPI) evaluated by
  a float filter, then expansion arithmetic; no constructions [snippet]. Lévy 2024 (arXiv 2405.12949)
  reports exact constructions up to 6x faster because they cache; also reports Cherchi 2022 crashing or
  wrong on some inputs [snippet].
- None of these predict cost before running. They all **escalate on demand** inside one decision. No
  published failure-rate numbers were found this session.

### Condition numbers and cheap ill-conditioning features
- Farouki and Rajan 1988 CAGD: condition of curve/surface intersections is the participants' conditions
  times factors from the normal directions; tangential intersections and singular points are inherently
  ill-conditioned [snippet].
- Bürgisser et al. 2017 SIAM J. Appl. Algebra Geom.: intersection condition number = 1/sin(alpha), alpha
  the minimum angle between tangent spaces; condition number theorem: 1/kappa = distance to the
  ill-posed set [snippet]. Bürgisser and Cucker, "Condition", Grundlehren 349 (2013), general framework
  [snippet, book unread].
- Sederberg and Meyers 1988 CAGD: normal-cone test; disjoint normal cones prove no closed intersection
  loops; cones come cheaply from control vectors [snippet].
- **Inference:** the angle feature measures how hard the problem is **in floating point**. In exact
  arithmetic a tangency is a zero sign, decided exactly. Catallax data agrees: P0 "line tangent" builds
  exact at 35,038 work and "near miss" at 3,625 (`checkpoint-1.md:45-60`) [read]. So 1/sin(alpha) is a
  good predictor of fitted-path cost and risk, and a poor predictor of exact cost on analytic pairs.

### Cost of exact arithmetic
- Yap, exact geometric computation: cost is set by algebraic degree and bit length of the expressions,
  with root separation bounds giving the worst-case precision [unread/memory].
- ESOLID (Keyser, Culver, Foskey, Krishnan, Manocha; SM 2002, CAD 2004): exact boundary evaluation
  restricted to algebraic degree <= 4 "for efficiency"; within one order of magnitude of a float system;
  uses lazy evaluation and float filters [snippet]. Ouchi and Keyser 2004: perturbations lengthen bit
  lengths and runtime grows sharply [snippet].
- Bicubic patch implicit degree 18; resultant degree of degree m and n surfaces is mn; coefficient bits
  grow roughly additively per resultant level [unread/memory].
- **Inference:** a theory-backed cheap exact-cost estimate is D * L^2 per step, where D = product of the
  participants' implicit degrees and L = limb count of the input coefficients. Plane 1, quadric 2,
  torus 4 keep D <= 8; fitted bicubics carried as 53-bit dyadic rationals give D = 324 per pair and
  large L. This separates Catallax's cheap and expensive cases better than "geometry kind" names do,
  and it is generic.

### Budgeted speculation with deterministic budgets
- Ski rental: optimal deterministic competitive ratio 2 (2 - 1/b discrete) [snippet]. Mapping: keep
  paying for the exact attempt ("rent") or switch to the fitted path ("buy"). Cap the attempt at the
  fitted cost F and total cost is at most 2x the better path.
- Purohit, Svitkina, Kumar NeurIPS 2018: a prediction plus a robust fallback; deterministic ski rental
  with hyperparameter lambda gives consistency 1 + lambda and robustness 1 + 1/lambda [snippet for
  framing, bounds from memory].
- Luby, Sinclair, Zuckerman 1993: universal restart schedules within O(log) of the best fixed cutoff
  [snippet for paper, unread/memory for bound].
- QSopt_ex (Applegate, Cook, Dash, Espinoza 2007): solve in double, verify exactly in rationals, raise
  precision (128, 192, 288 ... up to 3164 bits) on failure; each level gets a much lower iteration
  limit "because if any floating point precision fails, we can always move to the next" [snippet].
- Deterministic work limits in shipping solvers: Gurobi `WorkLimit` stops at the same point on the same
  hardware and settings; Z3 `rlimit` counts basic operations, adopted by F* and the Clang analyzer
  because timeouts are nondeterministic; caveat: rlimit tracks time badly on some inputs [snippet].
- **Inference:** the closest published practice to "try exact within N work units, else switch" is
  QSopt_ex's per-level iteration limit and Z3/Gurobi deterministic limits. Catallax's work units are
  the same idea, and `refund_on_refusal` is already the switch.

### Does any kernel choose its path per operation?
- Parasolid tolerant modelling (Jackson 1995): exact geometry by default; edges that cannot be made
  exact become tolerant (pcurves in a tube) [prior notes, read earlier]. That is per-result fallback,
  not a pre-run prediction.
- ACIS: TEDGE/TVERTEX created when exact intersection fails or on import; `tweak` re-intersects
  [prior notes].
- CGAL lazy exact: per predicate, on demand, not per operation [snippet].
- LEDA real / EGC: precision escalated per expression [unread/memory].
- ESOLID: fixed rule by degree (<= 4) [snippet].
- I found no kernel that predicts the exact-vs-float cost of a whole operation before running it.
  Practice is either a fixed class rule (ESOLID) or on-demand escalation/fallback (CGAL, Parasolid,
  QSopt_ex). [Inference from the sources above; absence not proven.]

## 2. Candidate designs

Notation: E = exact cost, F = fitted cost, both in work units.

| Design | Features read | Cost | Deterministic | Failure mode | Evidence |
|---|---|---|---|---|---|
| (a) kind rule | surface/curve kind of each participant | O(faces), negligible | yes | blind to cost inside a "cheap" kind (8.5 s STEP on cylinders/tori; 27-53M planar containment) and to cheap cases inside an "expensive" kind (k12 far-disjoint cut, 9M exact) | ESOLID degree rule; checkpoint-1 P0 |
| (b) feature estimate | implicit degree product D, input limbs L, bbox/normal-cone separation, 1/sin(alpha) at sampled edge points | O(control points); well under 1% | yes if features are exact or fixed-precision | thresholds carry theory but the constant does not; angle feature misleads on analytic pairs; misses uncharged costs | Farouki-Rajan, Bürgisser 2017, Sederberg-Meyers, ESOLID |
| (c) capped speculation | none before; measured work during attempt | the cap itself (wasted when it trips) | yes, if the meter is complete | wasted cap on every doomed attempt; expensive-but-finishing exact cases lost if cap too low; blind where the meter is blind | ski rental, QSopt_ex, Gurobi/Z3, Catallax `refund_on_refusal` |
| (c') memo per step | hash of step inputs | lookup | yes only with exact-input keys | exact-key memo adds little (unchanged steps are already checkpointed); fuzzy key changes results by history | Inference |
| (d) learned model | recorded work rows + features | small at run time | yes once frozen and versioned | ~40 rows from 2 parts: fits the test parts; path changes on retrain | Hutter 2014 needs large sets |
| (e) hybrid | (b) features choose the cap tier; (c) enforces it | (b) + cap | yes | combines (b) and (c) limits, but every misprediction is bounded by the cap | Purohit 2018 consistency/robustness |

## 3. Known cases

| Case | Measured | (a) kind rule | (b) D, L, bbox | (c) cap 20M, fallback if one exists | (e) |
|---|---|---|---|---|---|
| Box, 12 edges filleted | E 9,467-13,195; fitted network R1 72.3M (8 corners + 12 strips) | exact, right | D=1, exact, right | exact finishes, right | right |
| P1 bracket fillet-end vertex | E 13.6M whole build; F unknown (OCCT numeric; sample-part rolling edges 0.6-5.6M) | exact, 13.6M, under target | exact (assume analytic) | exact finishes under 20M; with a k x F cap and k=4, F~3M it would trip at 12M and then pay F: worse and loses exact | right with 20M tier |
| Loft later cuts k01, k04-k11 | E 81.6-98.9M final step, exhausted; no fitted Boolean exists | fitted: instant refusal; saves ~80-99M each | D large: same | trips at 20M, refuses; saves ~62-79M each (~560M over 8 cuts) | same as (c) |
| k03 rim notch | 24.3M then coverage refusal | saves 24M | saves 24M | saves 4.3M | |
| k12 far-disjoint | E 9.0M, builds, 22 faces | **wrong**: sends to a fitted Boolean that does not exist, loses the only success | bbox-disjoint promotes: right | right (9M < 20M) | right |
| M10-004..006 spline top cuts | 26-27M whole request, built exact | **wrong** (spline) | wrong unless tier allows | only right if the cap is lifted when no fallback exists | rule: cap binds only when a fallback exists |
| Owner rim fillet | exact builders refuse, refunded (cost unrecorded); rolling ball 6.0-6.4M | try exact, refuse, fall back: right | same | same; waste <= cap | same |
| Exact STEP export, curved faces | 807,063 work, 8.5 s; numerical STEP of rim 21.9-22.6M work | exact, 8.5 s | D=2-4, exact | cap never trips (Taylor arithmetic uncharged) | blind until `fit.rs` is charged |
| Analytic Booleans P0 | 2,813 to 3.17M; line tangent 35,038; near miss 3,625 | exact, right | exact if angle is NOT used to demote analytic pairs | right | right |
| Four-contact fan containment | B-rep 27.3-53.2M per query; field query 257,785 | exact B-rep (planes): costly | planes: same | trips at 20M, falls to field query | right |

Wasted work of wrong choices (Inference): choosing fitted on the box fillet costs ~72M extra and loses
exactness (worst single case). Choosing exact on the M8 loft cuts wastes ~80-99M each today. A 20M cap
bounds the second at 20M; nothing bounds the first except a rule or tier that says "try exact".

## 4. Recommendation

(e): theory-backed tier from (b), enforced by (c). Concretely:

1. Make the meter complete first. Charge `catallax-step/src/fit.rs` and every other uncharged rational
   path by limb product, as the decider already does. Until then no work-based predictor can see the
   STEP case.
2. Count refunded speculative work separately (`speculative_work` beside `used`) so machine cost of
   refused attempts is reported, not erased.
3. Wrap every exact attempt that has a fitted fallback in a sub-cap. Two tiers, from D * L^2:
   low (D <= 8 and small L: planes, quadrics, tori, closed-form blends) gets the full 20M target;
   high (fitted splines, dyadic 53-bit data) gets a small cap, e.g. 2M, enough for trivial cases.
   Cheap promoters: bbox or normal-cone separation moves a pair to the low tier.
4. Where no fallback exists, do not cap below the request allowance: a cap there only converts late
   successes (M10 top cuts) into refusals.
5. Use 1/sin(alpha) to set fitted-path density and tolerance, and to demote spline pairs only. Never to
   demote analytic pairs.
6. Skip the per-step memo (exact-key memo duplicates checkpoints; fuzzy memo breaks determinism) and
   skip the learned model until there are hundreds of held-out rows.

Why: no published kernel predicts whole-operation cost before running; the robust practice is
escalation under a deterministic limit (QSopt_ex, Z3 rlimit, Gurobi WorkLimit, CGAL filters). Ski
rental shows a cap gives a bounded worst case with no prediction at all; the prediction only tunes the
cap (Purohit). A wrong tier costs at most the cap, never correctness, because the fallback is the same
fitted path that runs today after a refusal.

### Done when
On a frozen corpus (every row of `docs/phase5-work-cases.md`, checkpoint-1 parts, box fillet, plus a
held-out ABC subset chosen before tuning):
- path choice and work counts identical at 1, 2, 4 and 24 threads and on two machines;
- predictor features cost <= 1% of the operation's work in every row;
- every case that builds exact today within 20M still builds exact (box, P0, P1, k12);
- refunded speculative work per operation <= its tier cap, and total machine work per row <= today's;
- wall time per million work units within a 10x band across all rows (catches uncharged arithmetic,
  e.g. STEP export);
- thresholds unchanged between the tuning set and the held-out set.

## 5. Not verified

- F for the P1 fillet-end vertex; geometry of that vertex.
- Cost of the refunded exact attempts in the rim cases (not recorded).
- Whether a fitted Boolean path will exist; without one the cut cases only gain faster refusals.
- Work units the STEP Taylor path would show once charged.
- The 95 to 448,000 units/ms spread (from prior notes).
- Pion's filter-failure numbers, Purohit's exact bounds, Luby's bound (snippet or memory only).
- The 2M and 20M tier caps are judgement, not derived; D <= 8 rests on ESOLID's degree-4 limit.
- No outside-CLI review: grok hit its usage limit.
