# Viability of spline-face rolling-ball fillets in the Catallax field (step 7, 4.2a)

Worktree Catallax repository at bb2f8375 (CX below).
Labels: [measured] = run this session; [read] = code or doc read this session;
[inferred] = my reasoning, not checked. Build output and logs:
viability-logs/
(offset_pair.log, spline_fillet.log). Machine: i7-9750H, 12 threads.

## 0. Verdict

Go with conditions. The spine half of the route works and is cheaper than the earlier study
assumed. The blocker is upstream: the field cannot hold a spline face at all, and its CSG of
implicit supports cannot hold one by the same route it holds analytic faces. Step 7 is
viable only after a new boundary-distance leaf for rational patches exists. That leaf is the
real cost and the real risk. Two corrections to the earlier route are also needed (section 7).

## 1. Prerequisite: certified distance to a spline face in the field

What the field has [read]:
- crates/catallax-field/src/import_spline.rs:1-21 admits only bilinear equal-weight patches,
  flat ones as planes, curved ones as a six-face "saddle roof" held by the exact quadric
  `gamma = alpha beta`. Everything else refuses `field.import-spline`.
- crates/catallax-field/src/sweep.rs:350-365: field lofts hold only straight translated
  correspondences and circular lofts; anything else refuses
  `field.loft-correspondence-unsupported` with "0 general loft patches admitted".
- The field builds solids as CSG of implicit supports. A bicubic patch has an implicit of
  degree up to 18, so the saddle-roof trick does not scale [inferred].
- Precedent for a different kind of leaf: crates/catallax-field/src/mesh.rs:1-11, a closed
  mesh held as "winding number one", sign from exact crossing count on a perturbed ray, value
  the true signed distance from the nearest triangle. 353 lines. A spline-boundary leaf would
  copy this shape: min over faces of certified point-to-patch distance, sign by certified
  ray crossings.

What exists to build on [read]:
- crates/catallax-shapes/src/distance.rs:1-5: deterministic interval branch-and-bound to
  compact trimmed analytic surfaces, with a caller trim predicate, `point_distance` (:676) and
  `face_distance` (:742). It refuses `Surface::BSpline` at :163 only because `Domain::new` maps
  a spatial box to parameter ranges through an analytic frame. For a patch the domain is its
  parameter rectangle already, and cell enclosures come from
  `BSplineSurface::evaluate_box` / `jet_box` (crates/catallax-shapes/src/bspline.rs:628, :645) or
  `RationalPatch::bezier_spans` / `enclosed_controls` (rational_patch.rs:788, :863).
- crates/catallax-shapes/src/spline_fillet/clearance.rs (348 lines): proves a point's closest
  point on a patch by a positive definite Hessian of |C - P|^2 on a box round the contact plus
  cell exclusion elsewhere. This is the local-minimum certificate a field distance needs. It is
  B-rep fillet code, so the field should rewrite it, not call it (section 5).
- crates/catallax-shapes/src/rational_patch/query/affine.rs (468 lines): line hits on a
  rational patch, usable for the ray-crossing sign.
- catallax-measure has spline bounds (src/bounds/spline.rs, finite_spline.rs), for prisms only
  [read headers].

Is it planned [read ROADMAP.md:196-198]: 4.2 "Field versions of fillets, sweeps and lofts",
done when "P5 builds in the field with no fillet refusal". P5 is analytic, so 4.2's done line
does not force a spline face into the field. 4.4 "Imported spline parts in the field", done
when "An imported spline part builds in the field". So general spline faces are 4.4 work, and
step 7 of 4.2a depends on a 4.4 deliverable. No 4.2a milestone doc exists in the worktree
(grep for "4.2a" finds nothing) [read].

Scope [inferred]: a whole-patch, untrimmed rational-patch boundary leaf (the same scope the
B-rep fillet accepts, M11-OFFSET-PROOF.md:184-200) is about 1.5-3k lines: patch distance
branch-and-bound with convexity certificate, ray-crossing sign, box lower bounds for display
and printing, typed refusals, controls. Trimmed faces (real imports) add certified
point-in-trim on parameter cells and roughly double it.

## 2. Measured: the M10 offset spines and the B-rep fillet

[measured] `cargo test --release -p catallax-shapes --test offset_pair -- --include-ignored`
and `--test spline_fillet`, CARGO_TARGET_DIR in the scratchpad.

| Edge | Spine work (M10) | Fillet work incl. spine | Spans | Bound |
| --- | ---: | ---: | ---: | ---: |
| M27 loft edge | 2,332,939 | 5,252,648 | 13 | 78 nm |
| Next loft edge | 2,666,401 | 6,970,934 | 16 | 5 nm |
| Biquadratic edge | 6,793,048 | 16,433,004 | 22 | 80 nm |

Identical at 1, 2, 4 workers and default; each spine is 1 open branch, 3 arcs, depth 3.
The 2.3M-6.8M claim is confirmed exactly. Between-caps variants: 8,199,428 and 9,685,010
units, 18 spans, 27 and 28 nm [measured].

Wall time [measured, release build]:
- `m11_offset_spines` alone: 2.77 s for 12 spine runs (3 edges x 4 worker settings),
  47.2M units, so about 17M units/s. Per edge by unit share: about 0.14 s, 0.16 s, 0.40 s
  [inferred split].
- `m11_fillet_faces` alone: 4.90 s for 12 fits, 114.6M units, about 23M units/s. Per edge
  about 0.22 s, 0.30 s, 0.70 s [inferred split].
- Whole spline_fillet test file: 14 tests in 26 s.

## 3. Does M10 accept analytic-against-spline offsets

No, not directly [read]. `OffsetPatch` wraps a `RationalPatch` only
(rational_patch/pair.rs:262-265), and `offset_meeting` (:287-299) passes `Source::Patch` for
both sides. `Source` is Patch, Atlas or Original (pair/cover.rs:8-12). No analytic variant.

Workaround [inferred]: the offset of a plane, cylinder, cone, sphere or torus is the same kind
of surface, so it is itself rational. Write the offset analytic face as a `RationalPatch`
(`RationalPatch::new(BSplineSurface)`, rational_patch.rs:238) with rational weights, e.g.
the circle `((1-t^2), 2t, 1+t^2)` whose Bernstein weights are 1, 1, 2, and pass offset 0
("A zero distance is the support itself", pair.rs:283). Limits: angular ranges need rational
half-angle tangents or an extended domain plus trim; a sphere patch must avoid poles; a full
revolution needs several patches. Not tested.

The B-rep's own spline-with-analytic fillets go through the numerical path
(spline_fillet/numerical_analytic.rs:1-3, "binary64 proposals ... carry no geometric
authority") [read]. So no certified spline-plane fillet exists to compare against.

## 4. Tangent and near-tangent supports

[read] M10-SPLINE-SPLINE.md:210-213: "An offset pair with a tangent contact refuses at the
split limit, as does a branch through a corner that both supports' borders share." The 24
M10 test pairs had 4 tangent contacts refused by name (ROADMAP.md:295). docs/research/tangent-
touch.md: the plane-patch tangent classifier was parked; border and curve touches undecided.

Where it bites in fillets [inferred]:
- Two offsets meet at the dihedral angle theta between the supports. Krawczyk contraction
  and the spine's Lipschitz constant both degrade like 1/sin(theta). theta -> 0 is a fillet
  that vanishes: an edge running out into a G1 join, a loft side meeting a cap tangentially,
  a variable-draft face. These are common in freeform parts (fillet chains propagate along
  tangent edges by default in commercial CAD), but I found no census in CX to count them.
- An exactly G1 edge has no fillet to build. Detect it and skip, not refuse.
- Multi-span supports put knot lines (internal G1/G2 joins) under the spine. M10 handles atlases
  (`Source::Atlas`), but the B-rep fit takes one-span patches only (M11-OFFSET-PROOF.md:184-200).

What would handle it [literature from memory, not re-read this session]:
- Krawczyk and interval Newton fail at singular roots by construction (singular Jacobian).
  Deflation restores a regular system: Leykin, Verschelde, Zhao, "Newton's method with
  deflation for isolated singularities", TCS 2006; Rump and Graillat, "Verified error bounds
  for multiple roots of systems of nonlinear equations", Numer. Algorithms 2010 (proves a
  double root of a nearby perturbed system); Mantzaflaris and Mourrain, deflation and
  certified isolation of singular zeros, ISSAC 2011.
- Tangential surface intersections: Ye and Maekawa, "Differential geometry of intersection
  curves of two surfaces", CAGD 16, 1999; Hu, Maekawa, Patrikalakis, Ye, CAD 29(9) 1997.
- Practical route that avoids the singular proof [inferred]: near a vanishing end the gap
  between the fillet and its supports is about r theta^2 / 8 measured along the bisector.
  Stop the certified spine at the station where sin(theta) is above a chosen floor, and in
  the tail answer the supports' distance with a bounded (not TrueDistance) value whose width
  is that gap. For r = 1 mm the gap is under 50 nm only when theta < about 6e-4 rad, so the
  tail is short, and the certified part must reach close to tangency. Whether the field's
  Krawczyk reaches that close is the open measurement.

## 5. Independence

Shared by the B-rep fillet [read]: M10 `offset_meeting` (topology: branch count, border ends),
spline_fillet's Krawczyk, jets, carrier, clearance. Generic primitives (`Ball`, `Decider`,
`BSplineSurface::evaluate_box` / `jet_box`, `RationalPatch`) are shared by every crate.

Recommendation [inferred]:
- Reuse primitives (catallax-numbers, evaluate_box, jet_box, RationalPatch, affine line hits).
  A bug there breaks many tests at once and is unlikely to hide.
- Do not reuse spline_fillet/*: write the field's own span certificate, ideally in different
  variables (e.g. contact on support A plus the closed-form or exact offset equation of B).
- The field does not need M10 at all for distance. A fillet's spine is the component that
  starts at the edge. The field can prove it by continuation: prove the two border roots with
  its own Krawczyk on the border system, then parametric Krawczyk span by span from one to the
  other (ParCont style). Completeness (no second branch) matters only for clearance, which the
  field proves anyway by its own point-to-patch distance. So M10 can be the B-rep's yardstick
  and the field stays independent.
- If M10 is shared anyway, a shared bug shows as both sides agreeing on a wrong branch. Catch
  it with a third check in tests: binary64 Newton on the offset system from a dense seed grid
  (e.g. 64 x 64 per support pair) plus brute-force sampled point-to-patch distances at 10^4
  points near the edge, compared to the field's enclosure. CX already uses this pattern
  ("independent binary64 Newton solve reproduces every end to six digits",
  M11-OFFSET-PROOF.md:43-45).

## 6. Cost estimates

Build per fillet [inferred from section 2]: the field's own spine by continuation should cost
less than M10 topology plus B-rep fit, because it needs no carrier, no Hermite rows and no
rail-domain proofs. Bound it by the measured B-rep: 5-16M units, 0.2-0.7 s per edge on this
machine. If M10 is used for topology, add 2.3-6.8M units, 0.14-0.40 s.

Per query [inferred, partly measured]: the field rolling_region panel cover (96 bits, at most
4096 panels, rolling_region.rs:14-15) on a closed-form centre curve: the 6 tests of
crates/catallax-field/tests/rolling_region.rs, which make at least about 100 classify calls
plus batched calls at 7 worker settings, finish in 0.04 s [measured]. So a closed-form
panel-cover membership query costs under 0.4 ms, likely tens of microseconds [inferred bound].
A spline spine needs per panel: a Taylor-model enclosure of C(tau) (cheap, stored), plus the
sector test against contact enclosures. The dominant cost is not the fillet but the
prerequisite leaf: each query also needs distance to the two spline supports and any nearby
face, a 2-parameter branch-and-bound per face. Expect 0.1-10 ms per query in exact/interval
arithmetic. Display and printing classify boxes, not points, so whole-part cost scales with
cell count.

## 7. Corrections to the earlier route

1. `rolling_region` is a membership classifier, not a distance. docs/M18-ROLLING-REGION.md
   [read]: "No exact distance, unique direction ... is implied by `Region::classify`". A
   `TrueDistance` answer needs the nearest point on the swept-ball surface, with the
   direction sector proved. Plugging into rolling_region gives sign, not distance. Decide
   which the consumer needs: printing and display need sign and box classification; mass
   and compare need distance bounds.
2. The earlier "eps <= min(100 nm, tol/2)" is reachable: the B-rep reaches 5-80 nm on the
   three edges with the same family of certificate [measured]. But the three edges are
   bilinear and biquadratic single-span patches (M11-OFFSET-PROOF.md:457-470). They say
   little about trimmed multi-span imports.

## 8. Risks

| Risk | Likelihood | Impact | Retire by |
| --- | --- | --- | --- |
| Field spline-boundary leaf is large and slow per query | High | Blocks step 7 | Spike: patch distance + ray sign on the M27 bilinear patch, time 10^4 queries |
| Near-tangent ends: Krawczyk stalls as theta -> 0 | Medium-high on real parts | Refusals at fillet ends | Sweep theta on a synthetic edge, record the smallest theta certified inside 20M |
| Analytic-spline pairs need a rational-patch wrapper | Medium | Spline-plane fillets refuse | Plane-vs-patch offset meeting through offset 0, one test |
| Shared M10 bug hides from both kernels | Low | Wrong fillet passes both checks | Field proves spine by own continuation; dense Newton third check |
| Trimmed and multi-span supports | High for imports | Scope limited to whole patches | Defer to 4.4; record as typed refusal |
| TrueDistance vs membership confusion | Medium | Wrong contract | Decide the value kind before coding |

## 9. Plan

1. Decide the value kind (membership plus bounded distance, or TrueDistance).
2. Spike (1-2 days of agent time): certified point-to-RationalPatch distance on the M27 patch
   pair, branch-and-bound plus Hessian certificate, written in catallax-field. Measure work and
   wall per query on 10^4 points. Retires the biggest risk.
3. Spline-boundary leaf for whole untrimmed patches (4.4 slice), ray-crossing sign.
4. Field spine by border roots plus parametric Krawczyk continuation, own code.
5. Swept-ball query on that spine, sector test, reach check.
6. Near-tangent sweep and tail policy.
7. Analytic-spline wrapper.

Done when: on the three M11 edges the field builds each fillet with its own proven
eps <= min(100 nm, tol/2) inside 20M work, identical at 1/2/4 workers; every B-rep display
vertex lies in the field enclosure widened by the B-rep bound; a dense binary64 third check
agrees; median query time recorded; a near-tangent control refuses or bounds by name.

## Not checked
- Per-query cost of a patch distance in CX (no existing code to run).
- Frequency of near-tangent fillet ends in real parts (no census found).
- Whether M10 accepts a rational patch written from an analytic offset (not run).
- The cited deflation papers were not re-read this session.
