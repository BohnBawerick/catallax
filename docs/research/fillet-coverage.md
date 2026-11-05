# Fillet coverage

This is the first coverage inventory, not M5 acceptance. The source baseline is
`86a81a5eaf35af18715978da066e6e8f209f1e75`. It includes the landed
[coaxial rounds](https://github.com/BohnBawerick/catallax/pull/298),
[collapsed-pole diagnostics](https://github.com/BohnBawerick/catallax/pull/299),
[fillets over fillets](https://github.com/BohnBawerick/catallax/pull/300) and
[equal-cylinder crossing construction](https://github.com/BohnBawerick/catallax/pull/304).
The older PLAN and M5 snapshots do not describe this baseline's full coverage.

## What counts

A build means the unchanged recipe CLI answered, wrote STEP, and the direct
OpenCASCADE import produced one raw-valid solid with positive finite volume.
The check does not sew, run ShapeFix, increase the import tolerance or retry a
smaller radius. It records maximum face/edge/vertex tolerance and rejects values
above 0.001 mm. That is an observer sanity limit, not permission to increase
Catallax's modelling tolerance. Matching OCCT volume within
`max(0.00001 mm^3, 0.0000001 * reference volume)` is a numerical cross-check,
not a distance certificate. The box control also has an analytic volume and
inside/outside witnesses.

`refused-with-reason` retains the actual code, detail and failing step. It is not
an impossible-geometry verdict. `wrong` means an invalid solid, failed material
witness, excessive tolerance or volume disagreement. A disagreement needs a
geometry/definition investigation before assigning blame to either kernel.
A crash, timeout, missing binding or missing response is `incomplete`, never a
clean refusal. OCCT's partial `BadShape` is not accepted as a result.

Three rows have no executable Catallax recipe. Two variable-radius cases use
a reference recipe plus an explicit OCCT endpoint law. The setback case retains
the attempted recipe, whose `overflow` field the current parser rejects, and
an independent station-plane construction. Their `schema.unavailable` records
come from the corpus, not from a kernel run; their clocks are `n/a`. Do not
count these as executed Catallax refusals or replace the laws with constant
rounds.

Construction, material side, naming/edit retention, whole-face deviation,
measurement, display, raw STEP and independent field evaluation are separate
gates. This first scoreboard checks the complete recipe/export route and raw
STEP, not every gate. `operations.json` preserves the construction spans even
when a later measurement/export refuses. No preview, interactive latency,
full-result distance, tangency or general self-intersection proof follows from
an OCCT topology check. In particular, `IsDone` plus `BRepCheck` is insufficient
for claiming arbitrary imported blends correct.

## Case inventory and public evidence

The table describes documented mechanisms, not measured parity. C names the
Catallax baseline; O, S and P refer to the sources below. A general surface
interface is not proof that every pair of surfaces succeeds.

| Family | Catallax baseline and executable controls | OpenCASCADE | SMLib | Parasolid |
| --- | --- | --- | --- | --- |
| Constant radius, convex and concave | Elementary plane strips, cylinder generators and circular reductions. `box-single`, `l-concave`, `boss-rim`, `hole-mouth`. | O1 constant-radius tangent contours. | S1 rolling ball with material side and selectable section. | P1 rolling-ball blend description. |
| Variable radius, including a zero end | No scalar-radius recipe spelling for a law. `variable-linear`, `variable-zero-end` explicitly expose the schema gap. | O1 endpoint, sampled parameter/radius and `Law_Function` overloads. | S1 variable rolling ball; documents zero radius at an end. | P1 documents varying radius along a blend. No local run. |
| Plane/plane at extreme material angles | Four polygon-prism controls approach 0 and 180 degrees without changing the requested angle or radius. Names say approximate degrees. | O1 numerical construction can fail; no universal angle guarantee. | S1 common solver/corner framework, no universal angle guarantee. | No public per-angle success boundary found. |
| Plane/cylinder, perpendicular and oblique | Circular torus and certified closed oblique bands. `cylinder-rim`, `oblique-cylinder`. | O1/O2 walking and approximation. | S1 trace rails, then create sections/surface. | P1 signed-offset centre path. |
| Plane/cone and cone/cylinder | Coaxial circular reductions landed. `cone-top`, `cone-base`, `narrow-cone`; `oblique-cone` is a noncircular gap. Cone/cylinder general contacts need a later control. | General fillet interface, not pairwise guaranteed success. | General surface solver, not pairwise guaranteed success. | General blend description, not pairwise guaranteed success. |
| Sphere neighbours | Coaxial circles versus varying-curvature contact. `sphere-dome`, `sphere-cylinder`. | O1/O2 general support/walking route. | S1 surface/topology-based construction. | P1 general rolling-ball description. |
| Torus neighbours and representation changes | Coaxial reductions include rational-radius plane/torus circles. `torus-slice` constructs before an export-budget refusal. `torus-intersection` keeps the unsupported irrational intersection-curve representation. | General support route; measure the actual case. | General support route; no per-case guarantee. | No public per-case guarantee found. |
| Spline neighbours and closed spline rims | Certified spline carriers do not imply general fillet construction. Not yet an executable row in this ordinary-recipe corpus. | O1 accepts general B-rep inputs; reported closed-spline failures in O4 are not locally reproduced here. | S1 NURBS support and approximate section controls. | General support claims, no local fixture or run. |
| Nonparallel cylinders, equal crossings and loops | Certified rings, saddles/loops and equal-crossing pole fans already exist. `cross-equal`, `cross-unequal`. | O2 walking can encounter domain, singularity and twist events. | S2 continuation, then corner dispatch. | Public internals insufficient to compare branch handling. |
| Chains, tangent chains and periodic seams | `box-mitre`, `tangent-boss`, `tangent-pocket`, and rotated boss seam. One selected edge must propagate along its intended tangent chain. | O1 `Add` constructs tangent contours. | S1 3x2 tangent continuation requires compatible radius/section. | P1 tangent-chain propagation. |
| Setbacks and preservation of neighbours | Version-4 defines station planes, but the parser does not yet accept `overflow`. `explicit-setbacks` records this schema gap, not a geometric refusal or automatic fallback. | No matching station-limit API used; the observer clips the removed region to the authored interval explicitly. | S1 separate radius, section and corner policies; not a claim of matching this recipe's station semantics. | P1 setback/collar controls are different authored policies, not silent retries. |
| Three or more selected edges at a vertex | Elementary equal planar spherical corners; mixed/curved corners remain bounded families. `box-all`, `l-mixed-all`; pyramid `vertex-4-one/two/all` covers valence four and partial selection. Valence five remains a corpus gap. | O1 warns that termination at 4+ incident edges is not handled. A successful local example does not remove that documented warning. | S1 NxN convex/mixed patch families; S2 dispatch is not complete for every N/M combination. | P1 vertex blends; no public general valence limit or tested parity. |
| Rolling over earlier fillets and unequal radii | Landed tangent cylinder/torus/sphere continuations. `unequal-ordered`, `mixed-ordered`, `roll-over-mitres`, `unequal-mitres`. | O1's limiting-face warning matters; do not use a failed OCCT round as a reason to reject Catallax's answer. | S1 tangent rollover and fillet/fillet trimming. | P2/P3 distinguish reblending and overflow capabilities from a universal promise. |
| Radius near face size, swallowed faces and mitres | Exact below/at/above controls retain roll-over semantics. `near-face-limit`, `at-face-limit`, `past-face-limit`, `swallowed-face`, `cylinder-collapse`. A regular-contact limit is not geometric impossibility. | O1 `StartSolFailure` can mean excessive radius for that construction; O4 reports threshold sensitivity. | S1 describes rail self-intersection, small adjacent edges, cliff rollover and feature interactions; S2 limits the shipped repair claim. | P1 smooth/cliff/notch outcomes differ. P2 specifically documents notch overflow on mixed geometry. |
| Partial rails and end walls | Closed oblique ellipse support is not support for clipped ellipse ends. `partial-oblique` preserves a real vertical end wall. | O1 warns when a blend escapes a limiting face. | S1 explicit end/corner handlers. | No public end-wall success guarantee found. |

## Sources and limits

- **O1:** OCCT 7.9.3 [MakeFillet header](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_3/src/BRepFilletAPI/BRepFilletAPI_MakeFillet.hxx),
  especially `Add`, `Build`, `BadShape` and `StripeStatus`. `SetContinuity`
  controls approximation continuity; it must not be advertised as an arbitrary
  G2 contact-to-support guarantee. `Build` explicitly lists 4+-edge termination
  and limiting-face escape restrictions.
- **O2:** [BRepBlend walking](https://occt3d.com/dev/doc/refman/html/class_b_rep_blend___walking.html)
  and the pinned [three-stripe corner builder](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_3/src/ChFi3d/ChFi3d_Builder_C3.cxx).
  Initial section, continuation, approximation and corner assembly are separate
  failure stages. No implementation code was copied.
- **O3:** [BRepCheck analyzer](https://occt3d.com/dev/doc/refman/html/class_b_rep_check___analyzer.html).
  Its optional exact curve-on-surface check applies only to SameParameter
  edges. The fast scoreboard does not enable it or claim that it proves
  self-intersection absence.
- **O4:** Upstream reports [1177](https://github.com/Open-Cascade-SAS/OCCT/issues/1177),
  [1371](https://github.com/Open-Cascade-SAS/OCCT/issues/1371) and
  [1695](https://github.com/Open-Cascade-SAS/OCCT/issues/1695). The independent
  reviewer supplied these as adversarial leads, not reproduced results. The
  last report is a reason to record tolerance growth, not merely `IsDone`.
- **S1:** NVIDIA's [SMLib filleting manual](https://docs.nvidia.com/smlib/manual/smlib/fillets/),
  sections Fillet Corners and Large Radius Filleting. It separates radius law,
  cross-section and corner generation. Topology-based G1/G2/G3 options must
  not be confused with the more limited surface/surface section implementation.
- **S2:** Pinned SMLib source `3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`,
  [corner dispatch](https://github.com/NVIDIA-Omniverse/usd-brep/blob/3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5/source/SMLib/src/SmFilletCorner.cpp#L1277-L1462)
  and [self-intersection handling](https://github.com/NVIDIA-Omniverse/usd-brep/blob/3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5/source/SMLib/src/SmFilletIntersector.cpp#L1631-L1718).
  The latter leaves its self-intersection test disabled and says the test is
  not implemented. Manual claims are not evidence of general automatic repair
  in this source revision. Repair can also change a circular section.
- **P1:** [Overview of Parasolid, v35 mirror](https://www.q-solid.com/Parasolid_Docs_V35/pdf/ov.pdf).
  Only indexed excerpts and the independent review were available; the full
  download timed out. Treat these capability descriptions as provisional, not
  source-verified solver internals or measured success. No fallback order or
  universal surface/angle coverage is inferred.
- **P2:** Siemens [Parasolid v34.1 release highlights](https://blogs.sw.siemens.com/plm-components/parasolid-v34-1-release-highlights/),
  section Blending Mixed Models. This is evidence for the named notch-overflow
  improvement, not every pure B-rep swallowed-face case.
- **P3:** Siemens [Parasolid v39.0](https://blogs.sw.siemens.com/plm-components/parasolid-v39-0-advancing-the-future-of-geometry-modeling/)
  documents changing a blend's size, shape and type during reblending. That
  does not demonstrate initial construction on our failing cases.
- **R1:** Farouki and Sverrisson, *Approximation of rolling-ball blends for
  free-form parametric surfaces*, 1996,
  [DOI 10.1016/0010-4485(96)00008-5](https://doi.org/10.1016/0010-4485(96)00008-5).
  Offset-intersection/contact-curve approximation suggests a construction
  family, not a guarantee of correct global trimming or a Catallax certificate.

The dot's baseline review was checked against the landed source and the SMLib
manual/source excerpts. Its upstream issue leads and Parasolid mirror notes
retain the qualifications above. Neither SMLib nor Parasolid has been run on
this corpus. There is no measured parity claim.

## Next constructions to compare

Keep [ownership round 1](../prototypes/fillet-ownership.md) and
[patch-first round 2](../prototypes/fillet-patch-first.md) parked. Round 1's
finite-sphere Boolean envelope spent its allowance without producing a body.
Round 2's fitted patches did not supply all production proof/trim obligations.
Neither result rules out contact-rail continuation with shared boundary curves.

Start with the small gaps that preserve current definitions: plane/torus
intersection-circle recognition, partial oblique end-wall clipping, and narrow
cone exchange. Compare their general contact/patch alternatives against the
same frozen recipes, not reduced replacements. Then tackle mixed corners,
unequal earlier mitres and swallowed-face continuation. The four-valence
pyramid must stay in the set so a special box-corner fix cannot masquerade as
general vertex support.

For each candidate, record initial-section solve, rail tracing/event handling,
corner assembly, trim ownership, certification and consumers separately. Reject
silent radius reduction, seam changes, healing, altered overflow policy or
unselected-edge rounding. Explicitly labelled alternative constructions can
be experiments, not replacements for the requested case. Keep losers here
with their failure stage and measured coverage. No new candidate is promoted
by this first PR. Investigate the oblique-volume disagreement below before
using either kernel's volume as a reference for another construction.

A fitted carrier must name its retained rolling-ball definition, local bound,
affected faces and cap exception. The existing policy targets 100 nm and caps
regular fitted fillets at 1 um, additionally bounded by neighbouring geometry.
The existing named P5 exceptions remain named exceptions. This work does not
raise those caps, relabel samples as whole-face proofs, or make approximate
answers `exact`. Any larger exception needs its own recorded reason and
review rather than an undocumented fallback.

## Measured baseline

The frozen 43-case run completed at 1,000,000 work units per executable recipe.
All status/code snapshots matched. Catallax delivered 13 checked exports,
returned 26 kernel refusals and had one volume disagreement. Three further
rows are schema-unavailable intents, not executed kernel requests. Thus the
table's 29 refusal records include three corpus records. OCCT independently
built 34 cases and refused nine. No process, binding or import was incomplete.
Fourteen Catallax exports imported as raw-valid solids; only thirteen met the
comparison checks. All thirteen accepted manifests report exact definitions.

The final selected fillet span answered in 25 cases, refused in 13, and was not
reached in five. Three of the latter are schema gaps and two are upstream
Boolean work refusals. These span counts are construction diagnostics, not
25 export, topology or whole-face proof passes. Export-budget failures must
not erase construction coverage already demonstrated by the logs.

The run used nice 10, affinity 0..11 and initial load averages
2.05/4.60/8.19. It took about 139 seconds including process startup and raw
imports, measured at the launching shell. It was not an idle-host benchmark.
Binary SHA-256:
`98696d0b1d262786c31fd217b3f26a9248b380e343b68b234155942fe12f7e9b`.
Original corpus SHA-256:
`8081b56b35eeeb154b62023c1b81e4d20769b0adaec921c50ca94aaf0001399a`.

### Stage accounting

The primary Catallax ledger separates construction from consumers. The raw
`counts` field remains unchanged for compatibility; `stage_counts` and each
row's `stage_result` locate the gap.

| Stage result | Cases | Meaning |
| --- | ---: | --- |
| checked-export | 13 | Raw STEP and numerical checks passed. |
| construction-refusal | 12 | The final requested fillet returned a non-budget refusal. |
| construction-budget | 1 | The final fillet exhausted its allowance. |
| upstream-budget | 2 | An input Boolean exhausted its allowance before the final fillet. |
| consumer-budget | 10 | The final fillet answered; later work exhausted the allowance. |
| consumer-refusal | 1 | The final fillet answered; STEP rejected its carrier. |
| failed-check | 1 | The oblique-volume disagreement remains unresolved. |
| schema-gap | 3 | The intent cannot yet be expressed as an executable recipe. |

Of the thirteen work-budget failures, ten are consumer failures, not missing
fillet construction. One is construction-limited and two are input-limited.
An earlier answered fillet in an ordered recipe does not qualify the final
requested fillet. An answered span is still not a topology, STEP or whole-face
proof pass.

The corpus now snapshots `expected_stage` as well as status/code. `--check`
detects a refusal moving between construction and consumers even if its
overall `numbers.operation-budget` code stays unchanged. Only snapshot
metadata changed; recipe definitions, OCCT selections and outcome targets did
not. The corpus with stage metadata has SHA-256
`9e0f92944a1b1c56b2c43a8f64a312355358996ba50f4cd35808d57ac380ab52`.

| Case | Catallax stage result | Catallax export | ms | Final fillet span | ms | OCCT | ms |
| --- | --- | --- | ---: | --- | ---: | --- | ---: |
| box-single | checked-export | built | 412.81 | answered | 0.64 | built | 7.40 |
| near-face-limit | checked-export | built | 480.30 | answered | 0.87 | built | 8.52 |
| at-face-limit | construction-refusal | refused-with-reason `fillet.unsupported-overflow` | 1.37 | refused | 0.28 | refused-with-reason `occt.fillet-not-done` | 7.77 |
| past-face-limit | construction-refusal | refused-with-reason `fillet.radius-too-large` | 1.49 | refused | 0.32 | refused-with-reason `occt.fillet-not-done` | 7.19 |
| swallowed-face | construction-refusal | refused-with-reason `fillet.radius-too-large` | 1.28 | refused | 0.26 | refused-with-reason `occt.fillet-not-done` | 8.21 |
| box-all | consumer-budget | refused-with-reason `numbers.operation-budget` | 41.68 | answered | 2.15 | built | 37.34 |
| box-mitre | checked-export | built | 4461.32 | answered | 1.83 | built | 15.80 |
| acute-5deg | consumer-budget | refused-with-reason `numbers.operation-budget` | 34.68 | answered | 4.11 | built | 8.03 |
| acute-halfdeg | consumer-budget | refused-with-reason `numbers.operation-budget` | 30.70 | answered | 3.62 | built | 7.12 |
| obtuse-174deg | checked-export | built | 28.94 | answered | 3.09 | refused-with-reason `occt.fillet-exception` | 2.31 |
| obtuse-179deg | checked-export | built | 28.04 | answered | 3.20 | refused-with-reason `occt.fillet-exception` | 2.28 |
| l-concave | checked-export | built | 1024.18 | answered | 0.65 | built | 19.12 |
| l-mixed-all | construction-refusal | refused-with-reason `fillet.unsupported-overflow` | 25.42 | refused | 1.76 | built | 87.67 |
| cylinder-rim | checked-export | built | 11452.36 | answered | 0.32 | built | 3.92 |
| cylinder-collapse | construction-refusal | refused-with-reason `fillet.radius-too-large` | 1.19 | refused | 0.14 | built | 4.15 |
| cylinder-past-collapse | construction-refusal | refused-with-reason `fillet.radius-too-large` | 1.12 | refused | 0.16 | refused-with-reason `occt.fillet-not-done` | 11.10 |
| boss-rim | checked-export | built | 8597.22 | answered | 0.30 | built | 11.32 |
| boss-seam-rotated | checked-export | built | 8583.24 | answered | 0.27 | built | 17.10 |
| hole-mouth | checked-export | built | 8468.40 | answered | 0.26 | built | 16.85 |
| cone-top | checked-export | built | 7.57 | answered | 0.65 | built | 3.73 |
| cone-base | checked-export | built | 8.08 | answered | 0.45 | built | 3.85 |
| sphere-dome | consumer-budget | refused-with-reason `numbers.operation-budget` | 4938.59 | answered | 0.38 | built | 21.90 |
| torus-slice | consumer-budget | refused-with-reason `numbers.operation-budget` | 49.56 | answered | 0.64 | built | 40.28 |
| torus-intersection | upstream-budget | refused-with-reason `numbers.operation-budget` | 24.19 | not-reached | n/a | refused-with-reason `occt.fillet-not-done` | 10.64 |
| cross-equal | consumer-budget | refused-with-reason `numbers.operation-budget` | 724.12 | answered | 96.78 | built | 96.42 |
| cross-unequal | construction-budget | refused-with-reason `numbers.operation-budget` | 1464.66 | refused | 1388.00 | built | 307.25 |
| oblique-cylinder | failed-check | wrong `volume.disagreement` | 3473.41 | answered | 3.10 | built | 33.85 |
| oblique-cone | construction-refusal | refused-with-reason `fillet.unsupported-contact` | 53.57 | refused | 0.06 | built | 60.90 |
| sphere-cylinder | construction-refusal | refused-with-reason `fillet.unsupported-contact` | 79.97 | refused | 0.10 | built | 212.03 |
| tangent-boss | consumer-budget | refused-with-reason `numbers.operation-budget` | 82.64 | answered | 4.37 | built | 21.91 |
| tangent-pocket | consumer-budget | refused-with-reason `numbers.operation-budget` | 1584.24 | answered | 5.10 | built | 31.24 |
| unequal-ordered | consumer-budget | refused-with-reason `numbers.operation-budget` | 97.54 | answered | 1.75 | built | 11.39 |
| roll-over-mitres | consumer-budget | refused-with-reason `numbers.operation-budget` | 35.76 | answered | 1.52 | refused-with-reason `occt.fillet-not-done` | 24.04 |
| mixed-ordered | checked-export | built | 4428.35 | answered | 1.81 | built | 24.89 |
| explicit-setbacks | schema-gap | refused-with-reason `schema.unavailable` | n/a | not-reached | n/a | built | 22.15 |
| variable-linear | schema-gap | refused-with-reason `schema.unavailable` | n/a | not-reached | n/a | built | 19.45 |
| variable-zero-end | schema-gap | refused-with-reason `schema.unavailable` | n/a | not-reached | n/a | built | 21.14 |
| vertex-4-one | construction-refusal | refused-with-reason `fillet.unsupported-end` | 23.99 | refused | 0.47 | built | 16.21 |
| vertex-4-two | construction-refusal | refused-with-reason `fillet.unsupported-corner` | 23.04 | refused | 1.97 | built | 108.26 |
| vertex-4-all | construction-refusal | refused-with-reason `fillet.unsupported-corner` | 28.37 | refused | 5.27 | built | 84.34 |
| partial-oblique | upstream-budget | refused-with-reason `numbers.operation-budget` | 95.93 | not-reached | n/a | built | 70.46 |
| unequal-mitres | construction-refusal | refused-with-reason `fillet.unsupported-corner` | 3.27 | refused | 0.13 | refused-with-reason `occt.fillet-not-done` | 17.54 |
| narrow-cone | consumer-refusal | refused-with-reason `step.unsupported-geometry` | 3.04 | answered | 0.38 | built | 3.29 |

### Separate higher-work probes

These four requests use 20,000,000 work units and a 45-second process guard,
not the default snapshot allowance. They separate coverage gaps from cheap-run
budget failures. They do not change the primary counts above.

| Case | Complete Catallax request | ms | Final fillet span |
| --- | --- | ---: | --- |
| box-all | built with raw-valid STEP | 13793.18 | answered |
| torus-intersection | `fillet.unsupported-contact` | 38.17 | refused |
| partial-oblique | `fillet.unsupported-end` | 137.14 | refused |
| narrow-cone | `step.unsupported-geometry` | 3.48 | answered |

The first two torus controls deliberately differ. The rational-radius slice
constructs and reaches export; the irrational circle represented by an
intersection curve reaches a contact refusal once its Boolean has enough
work. The narrow cone constructs its outer torus sheet but cannot export it.
These are not four instances of an impossible fillet.

### Unresolved volume disagreement

For `oblique-cylinder`, direct Catallax STEP import gives
618.494339781691 mm^3; OCCT's independently constructed fillet gives
618.511633890515 mm^3. Both are raw-valid. The difference is about
0.01729411 mm^3, above the 0.00006185 mm^3 comparison allowance. Keep this as
`wrong / volume.disagreement`, meaning a failed cross-check, not a confirmed
attribution to Catallax.

A diagnostic OCCT rebuild with
`SetParams(1e-6,1e-7,1e-7,1e-7,1e-7,1e-7)` reduced its recorded maximum
tolerance from 0.00005811 mm to about 0.0000001 mm, but its volume remained
618.511640323897 mm^3. Toggling knot-span integration and tightening numerical
volume integration from 1e-6 to 1e-10 did not close the gap. The existing
oblique regression's independent section value is near the Catallax result,
so OCCT is not an unquestioned volume oracle here. The default scoreboard
retains the default OCCT construction, without a silent tolerance retry.
Investigate the actual contact definition and exported geometry next.

The full runner intentionally exits nonzero for this disagreement even when
its expected snapshot matches. The six-test smoke/accounting suite passes;
that does not turn the failed geometric cross-check green.

## Rational spherical export frames

The export work after main `5e48ae71417b0b156a1f767eba6b736437140143`
targets the ten consumer-budget cases without changing their constructors or
the default file tolerance.

### Trials and decision

- Explicit file tolerances of 0.1, 1 and 100 nm produced no checked gains.
  An additional precision trial used 12, 8 and 6 significant digits at
  100 nm. Its preliminary per-decision policy differed from the CLI, so
  it is not production coverage evidence.
- Repeating the default 17-digit / 0.1 nm and 8- or 6-digit / 100 nm candidates
  with the production per-decision policy produced thirty request-budget
  refusals. No accuracy relaxation is promoted.
- Profiling actual CLI replays located the work. At 20M, `box-all`
  finished chart preparation at 6,990,760 work and export at 8,370,180.
  `sphere-dome` finished preparation at 525,143 and export at 2,508,825.
  `tangent-boss` still exhausted 20M during preparation. This is not one
  uniform surface-fitting bottleneck.
- Rational spherical frames produced a raw-valid `roll-over-mitres`
  export at 1M. Its volume differed from `3888 + 83*pi/3` by about
  0.000000000371 mm^3. The reported file bound was below
  0.000000004462 mm, with the default 0.0000001 mm request unchanged.
- Trying an existing clear-pole frame before pole-cap construction added no
  checked cases beyond the rational-frame result. Park that extra path.
  Keep the rational frames, with all previous diagonal candidates retained.

For perpendicular unit directions `a,r`, the direction
`3*a/5 + 4*r/5` is already unit. Rational unit source frames therefore need no new
square root before their chart and decimal proofs. The sphere's centre,
radius and point set do not change. The original axis permutations stay first;
the previous diagonals remain fallbacks for domains the new directions miss.

The new material regression checks the rounded block's top and upright
corners. Its footprint area is `396 + pi`; the first nine units have that
section. Integrating the final unit gives `324 + 56*pi/3`, hence
`V = 3888 + 83*pi/3`. This numerical volume control and five inside/outside
witnesses supplement raw topology checks; they are not a whole-part distance
certificate. The existing exact fillet definition is unchanged.

### Required certificate failures

A separate API probe used a 1M per-decision policy as well as a 1M request
allowance. For `sphere-dome`, the writer discarded a
`step.budget-exhausted` chart failure and returned an unseamed legacy file
that failed raw OCCT validity. The red regression reproduces that false
success without the fix. Local work, precision and decimal-certificate
failures now propagate instead of reaching that fallback; refusal leaves the
caller's output buffer unchanged. The ordinary CLI still builds this dome
with sufficient work. Do not count this safety correction as a new fillet
construction.

The updated expected export outcome and material witnesses have corpus
SHA-256 `61f27ce8dc3a0316a2bb5f83bf1676b3a21f1a9ab6bdedac7798fc76620f6175`. Recipe dimensions, selections and OCCT operations did not
change. Every kernel PR still reports a fresh full 43-case run at its final
committed head. CI's small live smoke is not that full run.

## Reproduce

The corpus is [parts/fillets/cases.jsonl](../../parts/fillets/cases.jsonl).
Each row stores a complete ordinary recipe or a clearly unexpressible law,
expected geometric intent, family tags, OCCT edge witnesses and the observed
baseline outcome. The OCCT adapter independently constructs stock, polygon
prisms and booleans; it does not import Catallax's unfilleted body. Edge
witnesses must select exactly one edge, never a face index or nearest guess.

Build the ordinary verification CLI, then use an external proof directory and
a Python environment with `cadquery-ocp==7.9.3.1.1`:

```sh
CARGO_TARGET_DIR=/tmp/catallax-fillet-target nice -n 10 cargo build --release -p catallax-cli --bins -j 2
nice -n 10 python3 -B tools/fillet-coverage/run.py \
  --binary /tmp/catallax-fillet-target/release/catallax-verify \
  --occt-python /tmp/occt/bin/python --out /tmp/fillet-coverage --check
CATALLAX_VERIFY_BINARY=/tmp/catallax-fillet-target/release/catallax-verify \
  /tmp/occt/bin/python -B -m unittest discover -s tools/step-regressions \
  -p test_fillet_coverage.py -v
```

`--case ID` repeats to select controls. The default allowance is 1,000,000
work units per recipe and the process guard is 30 seconds per request, not a
speed target. Each case freezes its input, command output, ordinary driver
bundle, independent reference and raw-import report. `scoreboard.json`
checkpoints after each case; `scoreboard.md` is a readable view. The report
records the reporting source's `source_head` and `source_dirty` flag separately
from the binary hash. Run the full corpus at each final PR head and quote that
head with its scoreboard; CI's smoke controls do not exercise all 43 cases.
Nonzero exit
means a wrong/incomplete outcome, or a changed baseline when `--check` is set.
A clean refusal is expected coverage debt, not a passing geometric target.
`--check` compares the recorded status/code and stage, not the timing; it does not waive
an expected `wrong` result. Use the same work allowance as the snapshot when
checking it. Require `status: completed` before treating a report as a full run.

The table's Catallax clock covers the complete CLI request. Its separate final
fillet span excludes later measurement/export. The OCCT clock includes its
construction and numerical checks, without Python startup. Raw STEP import
has a separate clock in JSON. These are single observations on a shared host,
not a paired speed comparison or application latency measurements.
