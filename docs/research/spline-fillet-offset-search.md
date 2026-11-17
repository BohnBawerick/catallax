# Spline fillets: the offset search, and the route that replaced it

Outcome: route changed. The direct search did not finish. A different route rounds one loft edge within 59 nm.
When: 11 to 12 November 2026.

## What we were trying to achieve

A fillet between two spline faces, on a loft and on an imported part, with the whole fillet face proved within 100 nm of its rolling-ball definition (hard limit 1 µm).

The centre of the rolling ball runs along the curve where the two faces, each moved inward by the radius, meet. Everything else in the fillet hangs off that curve, so the question was how to find and prove it.

## Why it was proposed

The spline-fillet milestone, now 2.1, was about to commit to fillets on arbitrary spline faces. The port scope study had said to scope free-form fillets separately from intersections. So a scope study was added on 11 November as a bounded experiment: does the cell search for spline cuts, which works for cuts, extend to the fillet's equations?

## What was tried

1. The direct extension: five equations in six unknowns, on one edge of a ruled loft.
2. A follow-up after the owner chose the 2026 engine: enclose the two normal lengths instead of searching them, which cuts the search to four unknowns. Tried on the same loft edge, the next edge and a curved control.

## Result

1. The six-unknown search refuses at proof depth 24 after 2,910,656 of 20,000,000 units. It does not run out of work. It runs out of depth.
2. The four-unknown search also refuses at depth 24 on all three edges. The loft edge uses 5,431,053 units.

No fillet face came out of either.

## Why the route was changed

The verdict in the note: this is evidence against that search, not proof that every exact offset search is too dear. It supports changing the proof route. It does not support raising the allowance, narrowing a tolerance or dropping small pieces.

The route that worked, from 12 November: give the spline-pair engine the exact offset faces as input, prove the centre curve there, then build the fillet face on the proven curve.

- The loft edge rounds at radius 0.1 mm in 10,349,665 units. The face is proved within 59 nm.
- It exports to STEP and opens as a valid solid in OpenCASCADE without repair.
- The same edge on a loft written by OpenCASCADE and imported rounds in 9,158,476 units, within 25 nm.
- Face and edge names survive a radius change.

## What is still open

The builder rounds one convex edge between two single-patch faces whose ends meet planes. Fillets across more than one patch, extra trims and junctions still refuse. That is milestone 2.1.

## Sources

- [M27-SPLINE-FILLET-SCOPE.md](../M27-SPLINE-FILLET-SCOPE.md): https://github.com/BohnBawerick/catallax/pull/400
- [M11-OFFSET-PROOF.md](../M11-OFFSET-PROOF.md): https://github.com/BohnBawerick/catallax/pull/445
- [M11-COVERAGE.md](../M11-COVERAGE.md): https://github.com/BohnBawerick/catallax/pull/521, https://github.com/BohnBawerick/catallax/pull/523, https://github.com/BohnBawerick/catallax/pull/526

## 16 November: numerical multi-patch body and raw STEP

The phase-2.1 continuation uses metered Newton offset contacts across atlas
charts, C1 cubic centre/section fits and paired world/native-UV end cuts.
The sampled tolerance is 1e-7 mm, not a whole-face or whole-body bound.
The perpendicular-plane seam and circular-section controls pass. An earlier
curved-loft carrier/cap-only probe passed at 1e-4 mm in 29,335,696 work; that
is separate from the current, stricter body result.

Declared numerical carriers use sampled chart-box coverage and nine UV-walk
locations per coedge. Only the new fillet carrier is declared numerical in
this probe. Source-carrier checks are unchanged. A shared edge incident to
a declared numerical face uses nine native parameters for its numerical
world-distance evidence, with the same 1e-7 mm tolerance. Exact-exact edges
retain their original clipped 3D/UV knot-span sample grid. The mixed-body
control rejects an exact-exact pcurve spike between the coarse locations;
the exact chart-gap and default native continuation/repair/support-edit/
handover controls also pass.

The smooth three-section loft's fillet body admits at 96,732,996 modeled
work within the 100M request allowance. The earlier 96,751,874 count also
included writing an unfilleted source STEP fixture. Independent circle-radius
and end-plane samples pass. Its measured final-admission work is:

| Check authority | Charged work |
|---|---:|
| Numerical face-local | 1,271,040 |
| Unchanged source-face-local | 6,540,517 |
| Cross-face checks and setup | 168,954 |
| Final-admission total | 7,980,511 |

These counts cover final admission only; the same body-build allowance also
covers the preceding source and fillet work. The STEP continuity fix below
lets the admitted body write and reread numerical STEP in a separate
80,265,942-work request. That is not a combined 100M end-to-end app result.
The kernel entry point is `catallax_features::fillet::numerical::build`;
ordinary recipe evaluation has not been migrated to this numerical result.

OpenCASCADE/OCP 8.0.1.1 reads the actual 644,190-byte STEP through
`StepToTopoDS_Builder`, with healing, `SameParameter` and `TransferRoots`
off. It is valid, with seven faces and one closed solid. A later box cut
is valid and reduces volume from 8.051606364758838 to 6.516691845842079 mm³.
The file SHA-256 is
`6f1984918902e4ee5f9ee2a57593aa33f9330f7b1c8f70ab52cef0c04fe3edaa`.

The [bounded attempt record](../phase5-work-cases.md#incomplete-21-bounded-attempts)
lists radius-edit, import, housing and all 13 frozen later-cut outcomes.
Only the far-disjoint M8 cut builds. No completed 2.1 milestone or stable
radius-edit names are claimed. The loft regression is
`m11_spline_fillet::numerical_fillet_builds_on_the_smooth_loft_across_its_station_seam`.

The interpolation and separate UV/3D-curve representation follow the
[SMLib rail-blend/lift workflow](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmFilletSolver.cpp#L8465-L8520)
and [OpenCASCADE's paired intersection curves](https://dev.opencascade.org/doc/refman/html/class_int_tools___face_face.html).
Numerical face-local and declared numerical-pair admission use Catallax's
finite sampling convention, not either kernel's complete wire certification.
No reference source was ported.

### STEP continuity encoding

Independent decimal rounding can split an exact C1 join whose continuity
is carried by pole relations at a degree-fold knot. The numerical STEP
helper now removes that one redundant pole strip before rounding, only
when the homogeneous inverse-knot-insertion equality is rational and exact.
The remaining knot multiplicity encodes C1 continuity structurally.
Non-removable joins are unchanged; no rounding tolerance is enlarged.

This follows the knot-multiplicity reduction used by
[OpenCASCADE `RemoveUKnot`/`RemoveVKnot`](https://occt3d.com/dev/doc/refman/html/class_geom___b_spline_surface.html),
restricted here to zero-error removal rather than tolerance-based smoothing.
[SMLib's `N_SrfRemoveKnotConditional` and boundary-constrained removal](https://docs.nvidia.com/smlib/manual/nlib/functions/index.html)
were checked as the alternative. No reference implementation was copied.
The decimal C1 regression and the existing `catallax-step` Rust suite pass.
Direct STEP from the admitted loft also passes its numerical reread and
the independent raw OpenCASCADE read/edit above. Optional diagnostic body
handover was removed after its separate 100M refusal; no transfer-engine
performance change or repeated-proof claim is made. Numerical body/history
layout 03 retains the new carrier declarations. Tagged v02 numerical body
and history records still read with empty declarations; unversioned
positional layouts continue to refuse.

### Review corrections

Output admission and endpoint projection use the stricter of the source
and requested fillet tolerances. A loose source therefore cannot admit a
new seam beyond the requested tolerance. This conservatively rechecks the
whole output at that stricter tolerance.

Native bindings are carried by face/boundary/use occurrence. Only replaced
rails, shortened border endpoints and inserted cap uses change; untouched
UV branches are retained, including two uses of a periodic seam edge.
Inserted cap bridges take their endpoint UVs from the two adjacent edited
occurrences, which were evaluated at their retained native edge parameters;
a separate projection of the new world vertex cannot replace those UVs.
Current-main surface kinds, continuous planar-circle clearance and
certified handover migration are preserved.


## 17 November: radius-edit end cuts

The radius change from 0.1 mm to 0.15 mm exhausted its 100M allowance in
the first end-plane cut. The diagnostic retained the original 1e-7 mm
tolerance and reported the typed geometric stage, not just the work limit.

[OCCT's `math_BracketedRoot`](https://occt3d.com/dev/doc/refman/html/classmath___bracketed_root.html)
uses derivative-free Brent interpolation inside a bracket.
[SMLib's fillet solver](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmFilletSolver.cpp#L2997-L3038)
instead sets bounded Newton equations for its coupled contact problem.
The end-cut candidate followed the OCCT bounded-interpolation approach,
using a secant proposal alternating with mandatory bisection rather than
the full Brent scheme. It still exhausted 100M in the first cut, so the
candidate was parked and the original bisection retained. No reference
implementation was copied.

In that experiment binary64 interpolation only proposed the next parameter.
Returned points still passed the same 128-bit plane-residual check, with
unchanged cap fitting tolerance and paired world/native-UV checks. The
measured request split is in the phase-5 record. No successful strict radius
edit or exact intersection authority was established.

