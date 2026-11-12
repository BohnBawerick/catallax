# P4 cylindrical probe proofs and allowance verdict

After merging main `64ad4de11ee379a15a96b632fa5eddf020888732`, the cylindrical
chart-proof change reduces P4's complete verification request from 1,720,668 to
**1,540,488 work units**, a saving of 180,180. Its manifest and STEP bytes are
unchanged. P4 still refuses at the 1M scoreboard allowance, during the rounding
boolean. The remaining complete-request gap is 540,488 units; M4 remains open.

The [initial profile](M4-P4-WORK-PROFILE.md) identified imprint probes and barrier
signs as the main cost. This attempt keeps their clearance proof while avoiding
some repeated expansion of world-coordinate equations. Tolerances and all three
work allowances are unchanged.

## The retained proofs

A cylinder chart interpolates its angular parameter `u(t)` and axial parameter
`h(t)` linearly, for `0 <= t <= 1`. Let `a` be its nonzero axis, `b` its radial
reference and `v = a cross b`. The chart denominator is
`b dot b + (v dot v) u(t)^2`, which is positive throughout the path.

Two boundary families therefore have simpler exact clearance proofs:

- A complete incident coaxial circle is the axial level `h = h_circle`.
  Excluding that linear level works for rational frames as well as algebraic
  frames. The old implementation restricted this proof to algebraic frames.
- A straight generator is an angular level. Its incidence and parallel direction
  are checked before using that level. Irrational generators may be stored as
  connected intersection curves rather than `Line` values. For a cylinder and
  a plane parallel to its axis, the component through the stored regular witness
  is one generator. The witness selects the branch; the other root of the plane
  section remains a different curve.

For a witness `w`, the finite angular level is

```text
u_generator = ((w - origin) dot v) / (b dot b + (w - origin) dot b)
              * (b dot b) / (v dot v)
```

A zero denominator identifies the generator omitted by this finite chart, which
the path cannot reach. Otherwise the existing polynomial exclusion proves that
`u(t) - u_generator` has no later zero. Source contact remains permitted only at
`t = 0`; a hit at the candidate endpoint still fails the test. A path contained
in a generator needs the existing finite-edge check. Unsuccessful recognition or
exclusion falls through to the existing world-equation method. Unresolved numeric
decisions retain the existing typed precision and work-refusal handling.

The boundary-arc direction and whole-path certificate still gate every accepted
probe. The change introduces no sampled clearance test and changes no surface,
meeting curve, trim or output tolerance. Existing ordered run-ahead execution
still handles independent barriers.

## Measured results

Both sides used the unchanged [canonical P4 recipe](../parts/P4-shelled-tray.json).
The baseline is `8c8d6628e55313514818d25659cfd32353630403`. The development
observer build used release optimization level 1 on both sides. These are work
counts, not a wall-clock qualification.

| Request phase | Baseline | With chart proofs |
| --- | ---: | ---: |
| Shell construction | 63,289 | 63,289 |
| Rounding boolean | 1,491,686 | 1,311,506 |
| Exact construction total | 1,554,975 | 1,374,795 |
| Result summary | 32 | 32 |
| STEP export | 178,485 | 178,485 |
| Complete verification request | 1,733,492 | 1,553,312 |
| Eager application build with 0.01 mm display | 1,819,042 | 1,638,862 |

Construction totals include the shell and boolean rows. The separate application
entry fits its existing 20M allowance; it does not qualify the unproved-preview
lifecycle. Verification at 20M and 100M gives the same answer. At 1M it refuses
with `numbers.operation-budget` after 999,999 units and leaves no STEP file.

The completed request's stage partition includes the shell's internal boolean:

| Stage | Work |
| --- | ---: |
| Intersect | 462,478 |
| Imprint | 726,001 |
| Classify | 182,173 |
| Sew | 309 |
| Other construction | 3,834 |
| Summary | 32 |
| STEP export | 178,485 |
| Total | 1,553,312 |

Removing the rational-frame restriction saves 95,082 units. Recognizing retained
generator branches saves another 85,098. A preceding power-to-Bernstein rewrite
saved zero units on P4 and was discarded.

At that isolated revision the remaining gap is **553,312 units**. Even exact construction alone exceeds
1M by 374,795. The attempted chart reductions do not close that gap. The current
coefficient enclosure caller still spends 205,620 units; circular-chart evaluation
spends 180,471 and root isolation 110,435. Further work must reduce additional
intersection and imprint proofs while retaining their certificates. The analytic
offset, shell and thicken paths remain the next separate capability landings.

## Validation after merging current main

The original candidate is `dcc640c8d005fdbd8b070160470eaf698a720fa5`.
Merging main `37c9920d660050dc80a5499050281cf3b6bd15f8` gives candidate
`57c7e2c886299baa750172022b40090baccb1991` with the same chart-proof logic.
The [fresh merge validation](evidence/m4-p4-main-merge.json) records these counts:

| Phase | Work |
| --- | ---: |
| Shell | 62,617 |
| Rounding boolean | 1,299,354 |
| Exact construction | 1,361,971 |
| Summary | 32 |
| STEP | 178,485 |
| Complete verification | 1,540,488 |
| Eager application build with display | 1,626,038 |

Main's other changes account for the additional 12,824-unit reduction. The
remaining complete-request gap is 540,488 units. At 1M the request still refuses
in the rounding boolean, after 999,995 units. No allowance changed.

All 70 Boolean library tests pass on that merged code. Its matrix used the
`work-span` observer, which compiles out kernel helpers. Those runs establish
observer determinism across settings. The production matrix below supplies the
thread and worker evidence. The STEP hash below is unchanged.

The rounded STEP comparison initially used the newer main baseline against the
older candidate branch. The isolated branch used 271,849 units without either
shortcut, 271,849 with generators only, and 271,801 with the circle shortcut or
both shortcuts. After merging main the candidate uses **271,513**, below the
paired check's main baseline of **271,561** at
`63a65a53733f8fa6dfd3d463ba6384960ad506e0`. The existing STEP comparison confirms
identical bytes, report fields, marks and deviation bound. No additional shortcut
change was needed. These are work counts, with no wall-clock claim.

## Verification and cases

[Portable evidence](evidence/m4-p4-probe-reuse.json) records source and binary
hashes, operation and stage counts, fresh request outcomes, thread comparisons
and the independent raw STEP check. Raw evidence retained privately.

All 68 Boolean library tests pass. The generator cases cover source contact, an
omitted chart pole, an interior crossing, a candidate endpoint hit, containment,
and a nonparallel line. They exercise translation, rotation and a non-unit axis.
Two intersection components sharing the same surfaces retain their distinct
witness-selected branches, including witnesses at another axial height. A
transverse circular section is rejected as a generator.

The initial request batches used the development observer. Their work and
geometry comparisons remain recorded in the portable evidence. Production
thread coverage is measured separately below.

The answered STEP SHA-256 remains
`f9a8b9a15a576a4989690f75e1bd9f697e052299b80b33c1bfb4fdc4634c7381`.
Direct `StepToTopoDS_Builder` import with pinned OpenCASCADE bindings, without
`TransferRoots` or shape repair, finds one valid solid with 15 faces. Its measured
volume remains 4,831.253636666568 mm³. This independent import measurement is not
a kernel volume certificate.

The protected 477-case replay is a separate unchanged CI gate. This artifact
does not claim a completed full replay, SolidDraft interaction, field comparison
or the complete M4 exit.


## Production validation on current main

[Current-main evidence](evidence/m4-p4-current-main.json) records the production
builds, with profiling features absent. The main snapshot is
`64ad4de11ee379a15a96b632fa5eddf020888732`; its merge into this branch is
`85a112faa66c4e73ed32fe947fb83eb5aa426287`.

| Complete verification | Main | With chart proofs |
| --- | ---: | ---: |
| P4 at 20M or 100M | 1,720,668 | 1,540,488 |
| P4 at 1M, refused | 1,000,000 | 999,995 |
| Near-miss control | 7,076 | 6,980 |
| Rounded boss native probe | 271,561 | 271,513 |

The answered and refused manifests match main. Strict P4 STEP bytes match main
and the hash recorded above. The unchanged rounded STEP comparator confirms
identical bytes, report fields and marks, with the same exact deviation bound.

All four requests agree at 1, 2, 4 and all 12 available kernel threads, using
1, 2, 4 and 1 persistent workers. Both normal and Boolean run-ahead modes pass:
32 fresh candidate requests, plus four main requests. Complete responses with
only clocks removed, including work and operation logs, agree across all eight
candidate configurations. Every worker exits zero. No cached result is used.

The old observer matrix did not exercise helper threads. It is superseded by
this production matrix. A separate production matrix also verifies all 14
analytic-shell cases at the same four thread and worker configurations.

The 1M gate remains a typed `numbers.operation-budget` refusal with no STEP
output. The complete request needs another 540,488-unit reduction. Exact
construction alone still exceeds 1M by 361,971 units. With the existing summary
and STEP costs, construction must fit 821,483 units; its current cost is
1,361,971. The rounding Boolean would need to fall from 1,299,354 to 758,866
units if the other phases stay unchanged. No further certified reduction is
claimed by this refresh. The preview lifecycle remains a separate open gate.


## Remaining-proof verdict

Fresh, separate development observers reproduce the production count exactly.
Their current partition includes the shell's internal Boolean:

| Stage | Work |
| --- | ---: |
| Intersect | 462,478 |
| Imprint | 713,177 |
| Classify | 182,173 |
| Sew | 309 |
| Other construction | 3,834 |
| Summary | 32 |
| STEP export | 178,485 |
| Total | 1,540,488 |

The eager application build with display remains 1,626,038 units. The profiling
feature disables helper execution; these stage measurements supply cost evidence,
while the production matrices above supply thread evidence.

Three retained proof callers account for 496,526 units: coefficient enclosures
use 205,620, circular-chart constant enclosures use 180,471, and root isolation
uses 110,435. Their combined cost is 43,962 units below the required saving.
With all other recorded charges fixed, reducing those three callers alone
cannot close the complete-request gap, even at a zero-cost lower bound. Other leading work includes 170,695 units in
nonzero geometric signs and 102,616 in exact-zero geometric signs.

The earlier power-to-Bernstein attempt saved zero units. The chart-level proofs
save 180,180, and no further certified reduction has been established in this
pass. The verdict is a failed 1M exit with retained proofs and a typed refusal;
it is not a claim that future exact proof reuse cannot reach the allowance.
Further work needs reductions across additional intersection and imprint proofs.
