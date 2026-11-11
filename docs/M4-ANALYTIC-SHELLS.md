# Exact analytic shell paths

This landing handles complete convex conical and cylindrical revolutions whose
planar caps were stored with reversed face sense. The earlier recognizers rejected
those caps before building their exact inset. Curved faces with reversed material
sense still require their concave envelope and retain a typed refusal.

## Distance proof and exact meetings

RECIPE sections 4.7 to 4.9 define the operations by metric distance and normal
sweeps. A parallel analytic carrier has zero deviation from its defining normal
offset when its metric parameters are changed exactly:

| Carrier | Exact change for signed distance `d` |
| --- | --- |
| Sphere | `radius += d` |
| Cylinder | `radius += d`, with unchanged axis |
| Cone with unit axis `a` and positive slope `s` | `apex -= d sqrt(1+s^2)/s * a`, with unchanged slope |
| Ring torus | `minor += d`, with unchanged major radius and core circle |

These identities apply to every point in the retained normal chart. They are
whole-carrier identities, not sampled deviation estimates. STEP decimal rounding
and pcurves retain the writer's separate certified error bound. Existing positivity,
nappe, collapse and spindle checks remain required. Finite faces, openings and
joins additionally need their projection and trim proofs. Translating infinite
supports alone does not establish a finite-face distance result.

The existing `CylinderSkin::cavity`, `ConeSkin::cavity` and
`CircularSkin::cavity` prove their finite domains and rebuild native exact meeting
curves. In particular, a conical opening is accepted only where its existing
finite kept-face projection proof applies. Opening its narrow end still refuses.
No sampled section, fitted intersection, tolerance increase or allowance increase
is introduced here.

A reversed planar face means that its outward normal is the negative of its
stored plane normal. Negating that normal produces the same plane and the same
outward halfspace with `Sense::Same`. Recognition uses this outward support while
retaining the original source loops and names. It does not reinterpret a reversed
curved face as convex.

For a closed shell, the proved cavity is an exact inner shell with reversed
material orientation. For a permitted planar opening, the source disk minus the
inset disk is the exact rim. The inset opening face is omitted, its reversed loop
becomes the rim's inner loop, and its native meeting edge is shared with the
reversed inset wall. Each opening must resolve to one source disk and one inset
disk on the same outward plane. Other topologies refuse with
`distance.analytic-rim-domain`.

Previously answered source domains retain their existing construction. The new
sewn assembly is selected for inputs whose reversed planar caps prevented the
previous recognizers from answering.

## Cases

The cases in `crates/catallax-features/tests/revolved_shell.rs` cover:

| Case | Test |
| --- | --- |
| Kept conical wall, floor corner, open rim and hole | `open_revolved_frustum_sews_true_inset_wall_to_the_closed_planar_rim` |
| Separate outer and inset cone apex vertices | `pointed_revolved_shell_keeps_distinct_outer_and_inset_apices` |
| Spherical cap and cylinder, with a finite nearest-face witness | `revolved_sphere_cylinder_shell_uses_the_finite_projection_proof` |
| Nonplanar opening requiring another assembly | `nonplanar_revolved_opening_keeps_a_typed_assembly_refusal` |
| Exact inward cone-cylinder meeting | `cone_cylinder_join_rebuilds_the_exact_inset_meeting` |
| Two open cylinder caps and their annular rims | `cylindrical_rims_share_the_inset_edges_at_both_openings` |
| Closed outer boundary and reversed cavity shell | `closed_revolved_shell_keeps_a_reversed_void_shell` |
| Consumed cavity and no kept faces | `consumed_cavity_and_no_kept_faces_follow_the_distance_definition` |
| Sphere offset, shell and front sweep with exact radii and material | `sphere_radius_paths_have_zero_analytic_deviation_and_the_defined_material` |
| Torus offset, shell and front sweep with unchanged core circle | `torus_radius_paths_preserve_the_core_and_have_zero_analytic_deviation` |
| Narrow conical opening without a finite projection proof | `opening_the_narrow_conical_end_keeps_the_finite_projection_refusal` |

The first case also runs after translation and rotation with a non-unit axis.
The tests check material on both sides of the inset, at the rebuilt meeting,
on the rim and inside its opening. They also check topology and the responding
operation path. Existing analytic tests retain the outward cylinder's toroidal
rim sectors, conical tip and cap cases, sphere normal-sweep folds and ring-torus
collapse restrictions.

## Measured capability

[Portable evidence](evidence/m4-analytic-shells.json) records the exact inputs,
source hashes, operation counts and artifact comparisons. Complete CLI requests
use the existing 20M allowance. No allowance or tolerance changed.

| Case | Shell operation | Complete CLI request |
| --- | ---: | ---: |
| Open revolved frustum | 5,654 | 41,166 |
| Cone-cylinder join | 9,435 | 42,021 |
| Open revolved cylinder | 108 | 3,934 |
| Posed frustum and non-unit axis | 6,319 | 49,899 |
| Closed frustum | 5,571 | 41,227 |
| Two open cylinder caps | 178 | 2,385 |
| Consumed cavity | 359 | 2,586 |
| Pointed cone and inset apex | 222 | 17,422 |
| Spherical cap over a cylinder | 668 | 1,531,469 |

The last case spends 1,528,886 units in STEP export after 2,583 construction
units. Its capability fits the existing allowance; this is no wall-clock claim.

All 53 selected feature tests and both torus STEP tests pass. Clippy passes with
warnings denied for the changed libraries and verification CLI. A consumer of the
production numeric library verifies nested path annotations, worker isolation,
disabled recording and zero accepted work for observations.

Four fresh 14-case CLI batches agree at 1, 2, 4 and all 12 available kernel
threads, using 1, 2, 4 and 4 persistent workers. Pointed and spherical-cap shells
also agree in separate four-thread-count checks. Work counts, manifests, logs
without clocks and strict STEP bytes agree. P4 and rounded-boss control answers
and STEP bytes match the retained baseline.

Direct `StepToTopoDS_Builder` import with pinned OpenCASCADE bindings, without
`TransferRoots` or repair, validates all 11 distinct checked solids. Independent
solid-decomposition volumes fit the unchanged reader-error rule. These numerical
import measurements are not kernel volume certificates. The unchanged protected
477-case replay remains a separate hosted CI gate.

## Full torus cavity export repair

The raw check exposed an existing full-torus shell export error. The retained
baseline produced the same invalid STEP bytes. A synthetic torus seam wire is in
the natural carrier sense; its bound must also follow the material face during
void-shell export. The old bound used shell reversal alone. The fix is limited
to loopless torus faces and is covered by a STEP regression test.

For major radius 5, outer minor radius 1 and inner minor radius 3/4, the defined
volume is `35 pi^2 / 8`. The old import reported an invalid solid with volume
154.2125687670208 mm³. The corrected file is valid and measures
43.179519254765765 mm³ against the independent 43.17951925476594 mm³ formula.
The error is also below the old numerical margin. Accepted work stays at 1,129.
No protected torus or revolved distance/subtraction inputs were found; the full
saved replay still runs in CI.

## Kernel-loop paths

An optional `path` field on an operation-log entry names the analytic
implementation that answered. Sphere and torus paths name their radius change;
cylinder, cone and circular-cell paths name their exact construction. New shell
paths end in `cylinder-sewn`, `cone-sewn` or `circular-sewn`. Recording asks no
numeric question and spends no accepted work. Nested operations keep their own
path, and disabled recording retains no annotation.

## Remaining domains

The independent field evaluator has its own source inventory and proof domain.
A narrow B-rep landing does not certify field parity or the full M4 gate.

Fresh baseline CLI probes at the existing 20M allowance found these additional
gaps before the analytic-shell change:

| Request | Work | Typed refusal |
| --- | ---: | --- |
| Inward offset of a sphere trimmed by a plane | 12,943 | `distance.unsupported-envelope` |
| Inward offset of a torus trimmed by a plane | 880,971 | `distance.unsupported-envelope` |
| Nonplanar opening of the revolved cylinder | 457 | `distance.analytic-rim-domain` |
| Narrow conical opening | 732 | `distance.cone-domain` |

The first two requests reached the offset operation and refused because their complete
curved trim and edge-envelope arrangement was unavailable. Full untrimmed radial
carriers already have their own exact paths. These trimmed cases remain separate
work, and spline-face offset remains outside the analytic proof domain.

Raw evidence retained privately.
