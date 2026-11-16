# M10 spline boolean qualification

Eleven original exchange requests now write checked approximate STEP: the
three side box/loft operations, through-all loft cut, M7 prism cut, three
imported-loft operations and three spline-pair operations. The three top
requests still refuse at 20M and produce no STEP file.

The request allowance is 20,000,000 work units, including construction and
export. Tolerances and existing proof checks are unchanged. A budget row is
the charged prefix before the next reservation fails, not the cost of a
completed proof. Larger or resumable proof jobs remain separate work.

## Shared fixes

An admitted positive-weight spline over a closed finite knot domain is
bounded. Its complete point set differs from an unbounded line, including
when every span lies on that line. The native spline comparison now answers
that case after its existing span checks. The separate overlap query still
refuses to replace an overlap interval with isolated points. The regression
uses an unregistered weighted, multi-span spline with a non-unit domain and
checks both operand orders. This argument applies to every admitted finite
native spline; it uses no part dimensions or fixture identities.

The recipe exchange now uses the existing manifest/5 spline vocabulary when
every spline face retains an exact `RationalPatch` or `RationalAtlas`
definition. Trimming and sewing preserve those definitions. Retained section
and patch-pair edges use the existing `intersection` tag. Fitted supports
keep their existing certificate routes. Summary vocabulary supplies no new
geometric authority: the writer must still prove its complete trims,
rounding, fitting and required pcurves before publishing STEP.

Previously the side cases, through cut and M7 cut stopped at summary with
`numbers.geometric-gap-unproven`. They now reach the writer and report its
actual missing proof or representation. Their B-rep construction work is
unchanged.

## Numbered refusals

Measurements compare main `37015127c32841e65419836c9d6b6c0c4d0577f9` with
kernel/CLI source `9d769ab42038779b283f1e51e5f098eed6294f63`.
The case generator is
[`test_m10_spline_boolean.py`](../tools/step-regressions/test_m10_spline_boolean.py).
Case numbers identify this inventory; native refusal codes are preserved.

| Number | Request | Before work | After work | Refusal after the change |
| --- | --- | ---: | ---: | --- |
| M10-001 | Side subtract | 8,175,141 | 8,219,519 | `step.unsupported-geometry` |
| M10-002 | Side union | 8,175,190 | 8,219,568 | `step.unsupported-geometry` |
| M10-003 | Side intersect | 8,171,004 | 8,214,023 | `step.unsupported-geometry` |
| M10-004 | Top subtract | 19,999,916 | 19,999,916 | `numbers.operation-budget` |
| M10-005 | Top union | 19,999,916 | 19,999,916 | `numbers.operation-budget` |
| M10-006 | Top intersect | 19,999,916 | 19,999,916 | `numbers.operation-budget` |
| M10-007 | Through-all loft cut | 2,526,697 | 2,527,978 | `step.unsupported-geometry` |
| M10-008 | M7 quadratic prism cut | 1,802,159 | 1,802,159 | `numbers.geometric-gap-unproven` |
| M10-009 | Imported loft subtract | 115,984 | 115,984 | `intersect.parametric-coincidence` |
| M10-010 | Imported loft union | 115,984 | 115,984 | `intersect.parametric-coincidence` |
| M10-011 | Imported loft intersect | 115,984 | 115,984 | `intersect.parametric-coincidence` |
| M10-012 | Two spline lofts subtract | 74,693 | 74,693 | `intersect.unsupported-pair` |
| M10-013 | Two spline lofts union | 74,693 | 74,693 | `intersect.unsupported-pair` |
| M10-014 | Two spline lofts intersect | 74,693 | 74,693 | `intersect.unsupported-pair` |

M10-001 through M10-003 spend 44,378, 44,378 and 43,019 more units on
summary/export after successful construction. M10-007 spends 1,281 more.
These are the slower requests in this matrix. They stop at the writer's
`geometry has no specified STEP form` boundary. M10-008 reaches the complete
rational-patch check, which requires one boundary loop.

The imported requests supply the existing independent OpenCASCADE fixture
`m27-loft-occt.step`, with its SHA-256 checked against the actual bytes.
It has four bilinear spline side faces and two planar caps. Import completes
in 114,370 units; the Boolean then needs a finite curve-equivalence proof.
The two-spline requests build both operands, then stop at the atlas/atlas
surface-pair coverage boundary. These early refusals establish neither
classification nor sewing of those requests.

## Numerical STEP embedding

The exact writer has no STEP form for a multi-chart atlas face or a
patch-pair curve, so it refuses `step.unsupported-geometry` for all seven
constructed requests. Their Boolean results keep complete finite trims, so
this is a representation gap only. Following the ROADMAP "Numerical first"
rule, such a refusal now tries a numerical embedding at the exact writer's
1e-7 mm curve tolerance. Each atlas face is written as its retained source
surface. Patch-section and patch-pair edges become fitted cubics, and every
edge use on a spline face gets a fitted native UV curve. Numerical admission
and the writer's decimal payload are checked at that tolerance. Answers are
marked approximate, with no complete body bound. See
[the interface note](interface-changes/2026-11-16-numerical-step-embedding.md).

Work-PC replays at 20M, one worker:

| Number | Request | Total work | Embedding and STEP | Raw OCCT faces/edges/vertices | Volume mm3 |
| --- | --- | ---: | ---: | --- | ---: |
| M10-001 | Side subtract | 9,264,694 | 1,015,408 | 9/18/12 | 10.700529128388986 |
| M10-002 | Side union | 9,266,165 | 1,016,830 | 9/18/12 | 30.700529128388986 |
| M10-003 | Side intersect | 8,943,571 | 699,781 | 6/12/8 | 3.299479539097737 |
| M10-007 | Through-all loft cut | 3,363,176 | 819,153 | 8/18/12 | 13.00000866748672 |
| M10-012 | Two spline lofts subtract | 4,146,852 | 556,964 | 5/9/6 | 5.140628182592781 |
| M10-013 | Two spline lofts union | 4,178,287 | 587,495 | 5/9/6 | 19.1406368500795 |
| M10-014 | Two spline lofts intersect | 5,228,555 | 547,923 | 5/9/6 | 8.85938048489394 |

Every file opens in OpenCASCADE (OCP 8.0.1.1) through the repository's
healing-off reader as one valid solid with the manifest's counts, and a
later box cut stays valid. The volumes agree with each other: subtract plus
intersect gives the loft volume 14.0000087 in both families, and union minus
subtract gives the side box volume of 20.

## Verification

The finite-query solver moves M10-009 through M10-011 past their native
curve/patch gate. Positive-weight control hulls prove finite exclusion. An
invertible constant-weight bilinear patch admits a quadric pullback for
native quadratic curves and lines. Roots keep their original parameters;
multiplicity and the two finite germs distinguish tangency from crossing.
Whole finite containment has a separate control-hull check. Partial overlap
and ruling lines still require explicit interval ownership and refuse.

The exact writer's finite patch-row gate remains an exact refusal. Its
`numbers.geometric-gap-unproven` now tries the shared numerical embedding,
with the exact UV-domain preservation from the embedding owner. The result
stays numerical at 1e-7 mm, with source and emitted checks and no complete
body bound. The existing publication/reread tail is a shared helper; this
does not change its file checks or the `embedded_exact` signature.

The end-to-end test runs fresh requests with request-worker/kernel-thread
counts 1/1 and 4/4. It compares manifests, operation logs and work, excluding
diagnostic timings and statistics. Eleven cases must answer with approximate
STEP after their retained exact-writer refusal. The three original top
requests keep their classification budget refusal. A replay of M10-008 over
a seeded stale file must replace it with the same checked STEP bytes.

### Additional original-request acceptance

At source `f6944020c61064edfd816efb18a36620d46f1226`, each unchanged original
`catallax-step/1` request completes construction, numerical admission, STEP
and disk reread within the unchanged 20M allowance. The four actual files
pass the healing-off OCP 8.0.1.1 reader, topology and material checks, volume
comparison and a later cut. The M7 physical hole remains present. See the
[generated receipts](phase1/spline-boolean-finite-queries/acceptance.json).
The per-body importer adds 353 work units to each imported request; all four
STEP files remain byte-identical to the raw-checked predecessor files.

| Number | Total work | Faces/edges/vertices | Volume mm3 |
| --- | ---: | --- | ---: |
| M10-008 | 2,240,515 | 9/18/12 | 14.333333333333336 |
| M10-009 | 4,461,368 | 10/24/16 | 6.695312500000001 |
| M10-010 | 4,695,911 | 12/30/20 | 22.6953125 |
| M10-011 | 3,829,629 | 6/12/8 | 0.9713541666666666 |

### Numerical top source construction

The separate numerical API reuses complete exact intersection and imprint,
then classifies material witnesses against planar stock or finite rays on
the retained loft charts. Every reported ray hit keeps its native UV,
normal enclosure and original-source residual checks. Numerical labels
never enter exact classification or incidence facts. Selected face uses
form a private surface collection with one opposed pair per edge; only an
admitted `NumericalBody` leaves the embedding callback.

The focused M10-006 source test independently checks all fifteen witnesses
against the recorded natural-cubic loft and stock laws, then admits the
four-face, six-edge, four-vertex numerical body at 1e-7 mm under 20M.
This is source-body qualification, not a completed STEP request. Direct
M10-004 and M10-005 numerical routes exhaust during native-trim admission.
M10-006 admits its source at 19,908,212 work, then exhausts during STEP at
19,999,962. The original top CLI path is unchanged; the unqualified slower
adapter remains private. No top STEP or whole-milestone completion is claimed.

The method follows OCCT's finite curve/quadric route in
[`IntCurveSurface_Inter.gxx`](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/IntCurveSurface/IntCurveSurface_Inter.gxx)
and nearest regular finite-ray transition in
[`BRepClass3d_SClassifier.cxx`](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/BRepClass3d/BRepClass3d_SClassifier.cxx).
The [SMLib port map](M25-INTERSECTION-PORT-SCOPE.md) places the bounded
curve/surface query and finite overlap classification in step 6. No OCCT
code was copied.

The feature test retains the side Solid, one-material-lump and structural
checks. Its through-all cut also requires one solid and two cap openings.
The all-solid top target remains ignored and unfulfilled.

```sh
cargo test --release -p catallax-intersect --lib spline_line::tests
cargo test --release -p catallax-features --test m10_box_loft --test later_cuts
python3 -B -m unittest discover -s tools/step-regressions -p test_m10_spline_boolean.py -v
```

The existing Boolean benchmark covers 49 cases, two poses and three
operations: 294 whole runs and 1,764 stage rows. Before and after, every
recorded decision/work row agrees. Whole-run decider work totals
1,399,726,972 in each run. No row gets slower in deterministic work. These
are benchmark decider counters, not complete CLI request counters or a
wall-clock latency claim.

The thirteen existing `p5-later-cuts` recipes also ran as fresh seamed STEP
requests at 20M each. Comparable manifests, strict STEP bytes and every
request's work agree before and after. Total charged request work is
213,818,257 in each run, including refused prefixes. No cut gets slower.

Raw commands, binary identities, request/response files and benchmark logs
are retained outside the repository. The 477-case saved-answer replay remains
a required forge check; this capability matrix does not substitute for it.
