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

## 16 November: numerical multi-patch body, STEP still refused

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

The smooth three-section loft's fillet body admits at 96,751,874 work within
the 100M request allowance. Independent circle-radius and end-plane samples
then pass. Its completed final-admission work is:

| Check authority | Charged work |
|---|---:|
| Numerical face-local | 1,271,040 |
| Unchanged source-face-local | 6,540,517 |
| Cross-face checks and setup | 168,954 |
| Final-admission total | 7,980,511 |

These counts cover final admission only; the same request allowance also
covers the preceding source and fillet work. Numerical STEP construction
then refuses with `step.numerical-rounding-topology`: “the written atlas
source no longer has the original charts”. No filleted STEP, radius-edit
result, imported-part fillet or housing acceptance is claimed. The test is
`m11_spline_fillet::numerical_fillet_builds_on_the_smooth_loft_across_its_station_seam`.

The interpolation and separate UV/3D-curve representation follow the
[SMLib rail-blend/lift workflow](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmFilletSolver.cpp#L8465-L8520)
and [OpenCASCADE's paired intersection curves](https://dev.opencascade.org/doc/refman/html/class_int_tools___face_face.html).
Numerical face-local and declared numerical-pair admission use Catallax's
finite sampling convention, not either kernel's complete wire certification.
No reference source was ported.

