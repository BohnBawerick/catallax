# M10 spline section imprint queries

This slice connects certified section branches to Boolean curve meetings and
finite atlas incidence queries. It does not yet prove the box-with-loft
Boolean, complete trimmed-atlas classification or sew its result. Section
branches still have no STEP writer.

## Exact curve meetings

A section branch retains its cutting surface, finite charts and complete
component cover. Every meeting with another curve lies among that curve's
complete meetings with the cutting surface. The Boolean keeps a candidate
only after proving membership in the retained section component.

For a finite B-spline curve, each knot span supplies homogeneous polynomials
`X, Y, Z, W`. Substituting the cutter's implicit equation of degree `d` gives
`W^d F(X/W, Y/W, Z/W)`. Positive weights clear the denominator without losing
a root. Rational coefficients use the existing complete univariate root
solver; every retained root belongs to its closed knot span. Shared-knot
meetings are deduplicated exactly. Other coefficients refuse by name.

When the other curve satisfies the cutting equation throughout, a second
analytic constraint can supply a complete finite candidate set. Each
candidate must still lie on that curve exactly. A missing finite constraint,
an unresolved overlap or containment over only some spline spans returns
`boolean.patch-section-meetings`. Candidate incidence must be exact.

The six focused regressions cover:

- A plane section of the exact loft-like patch against its cylinder section,
  in both input orders. They meet only at `(0, 0, 5/2)`. Moving one patch by
  10 along y leaves no meeting, even though the cutting coordinates agree.
- A line lying in the section's cutting plane, with both exact end meetings.
- Reversal of the same certified component.
- A spline meeting at a shared knot, retained once, and typed refusals for
  whole-span and partial-span containment without a complete finite query.
- Exact finite atlas incidence, including public face-membership exclusion
  for a point inside the control hull but off the carrier.

## Finite incidence before body queries

An atlas has no globally defined implicit side. Incidence instead asks every
finite chart for an exact preimage. Public face membership can therefore
exclude an off-carrier point without asking for a side. Body queries on
finite rational patches and atlases use incidence before the existing ray
parity solver. This does not supply the missing trimmed-atlas solver.

One later body query at `(2, 1, 1)` reaches
`intersect.atlas-section-row`: its section contains a whole chart boundary
row. That refusal remains open; this slice does not turn it into membership.

## The loft Boolean still refuses

The recipes and complete section checks are in
[`m10_box_loft.rs`](../crates/catallax-features/tests/m10_box_loft.rs).
Subtract, union and intersect have the same first refusal for each box.
The preserved main is commit `53da1a599cf3b50ded26bdb1cf3eca4d49491603`.
Fresh 1-, 2- and 4-worker CLI runs agree on outcomes, details and work.
The application allowance remains 20,000,000 units.

| Box | Preserved-main refusal and work | This slice's refusal and work |
| --- | --- | --- |
| Side: 2 x 10 x 1 at (1, 1, 1/2) | `boolean.imprint-ambiguous`, unsupported elimination chart; 1,074,776 | `shapes.parametric-query`, cross-kind spline support equality; 9,488,794 |
| Top: 10 x 10 x 1 at (-3, 5/2, 1/2) | `shapes.parametric-query`, missing trimmed-body side query; 2,449,104 | `numbers.operation-budget`; 19,880,186, exhausted |

These are advances in the queries consumed by imprinting, not proven solids.
Cross-kind support identity and the closed-branch path still need work under
the same allowance. Classification, sewing and section-branch STEP output
remain required for M10.

## M8 cuts and STEP readback

The [M8 list](M8-LATER-CUTS.md) was rerun at the unchanged 20M application
allowance. The disjoint-tool cut builds in 12,932,203 units. The other twelve
cuts retain named missing-intersection or budget refusals. No additional M8
cut builds in this slice. All thirteen cuts and canonical P5 have identical
outcomes, details, work and STEP bytes at 1, 2 and 4 workers.

| Case | Outcome | Missing step | Units |
| --- | --- | --- | ---: |
| k01 wall hole | `numbers.operation-budget` | Exact query within the allowance | 19,999,958 |
| k02 top slot | `numbers.operation-budget` | Exact query within the allowance | 19,997,792 |
| k03 rim notch | `intersect.parametric-surface` | Circle against fitted rolling-envelope spline | 2,857,959 |
| k04 torus pocket | `intersect.traced-pair-required` | Skew cylinder/torus tracing | 3,295,244 |
| k05 plane notch | `intersect.unsupported-pair` | Fitted rolling-envelope spline with plane | 3,438,651 |
| k06 cylinder drill | `intersect.unsupported-pair` | Fitted rolling-envelope spline with cylinder | 3,439,512 |
| k07 cone drill | `intersect.traced-pair-required` | Non-coaxial cone/cylinder tracing | 3,934,213 |
| k08 sphere dimple | `intersect.traced-pair-required` | Non-coaxial cone/sphere tracing | 2,835,452 |
| k09 ring groove | `intersect.unsupported-pair` | Fitted rolling-envelope spline with torus | 2,837,814 |
| k10 axis section | `intersect.unsupported-pair` | Fitted rolling-envelope spline with plane | 2,832,545 |
| k11 offset section | `intersect.unsupported-pair` | Fitted rolling-envelope spline with plane | 3,433,193 |
| k12 disjoint tool | Solid | None | 12,932,203 |
| k13 fillet first | `numbers.operation-budget` | Exact query within the allowance | 19,997,484 |

The disjoint cut and canonical P5 both export STEP. The independent
OpenCASCADE reader opens each as a valid shape with 22 faces, matching the
model. This is no claim about STEP output for the refused loft Booleans.

## Saved-answer replay

Fresh preserved-main and candidate replays agree on all 477 saved cases:
normalized manifests, exported bytes and work records are unchanged. Both
runs have 457 answers, 8 `boolean.self-touching-shell` refusals,
12 `numbers.operation-budget` refusals and no crashes. The saved corpus
retains its existing 100M request allowance. The loft and M8 capability
checks above use 20M.

Raw evidence retained privately.
