# M10 spline boolean qualification

The three side box/loft operations and the through-all loft cut construct
solids. All fourteen exchange requests below still refuse. None produces a
STEP file, so there is no new OpenCASCADE export-validity claim.

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

## Verification

The end-to-end test runs fresh requests with request-worker/kernel-thread
counts 1/1 and 4/4. It compares complete manifests, operation logs and work,
excluding diagnostic timings and statistics. It requires the five constructed
cases to reach STEP, import to complete for M10-009 through M10-011, and every
refused output directory to contain no `answer.step`. A replay over a seeded
stale STEP file must remove that file and retain the same outcome and work.

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
