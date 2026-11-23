# Fitted sweep comparison

The fitted compare adapter uses the application's numerical sweep
admission. It follows the existing SMLib sample-error-driven spline fit
and OpenCASCADE's point-pair guide direction under the saved
equal-native-parameter law. The independent field evaluates the authored
motion with outward intervals. The two routes share validated profile
charts and arithmetic. [SMLib fitting](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmCurve.cpp)
· [OpenCASCADE pipe shell](https://dev.opencascade.org/doc/refman/html/class_b_rep_offset_a_p_i___make_pipe_shell.html)

For a line perpendicular to the sketch, the carried section planes are
parallel and distinct. A projected line guide must stay separated from
zero, as proved at the exact minimum of its squared length. For a
rotation-minimising quarter helix the field uses Bishop transport. Its
radial and binormal axes remove the helix's integrated torsion. A
positive separation test compares tangent advance against sketch tilt
and the maximum normal rotation at the profile's radius. It also
certifies transverse motion at each section. An authored twist rotates
about the transported tangent. [Bishop's parallel frame](https://doi.org/10.1080/00029890.1975.11993807)

The comparison associates each fitted side with its source segment and
native parameter rectangle. It checks that the fitted spans cover the
source interval once per chart. On each span it bounds the difference
of corresponding longitudinal control curves. The linear-interpolation
remainder on a unit interval is at most one eighth of the maximum
second derivative. That term plus enclosed endpoint errors bounds a
span. The solver subdivides dyadically when the bound is too wide.
Positive rational weights turn the rail bounds into a bound on the
whole section. A separate relative-weight term accounts for rounding
of the fitted chart weights.

At each end the shared trim rails bound the difference of the cap
outlines. Projection onto the authored cap plane preserves the projected
cap's winding away from the swept outline. The fitted height over that
plane attains its extremes on the cap boundary. Combining the planar
and height bounds gives the cap bound. The display comparison adds the
certified chord bound for each face and uses the maximum of those sums.

The adapter covers perpendicular line spines with optional line guides
and twist, and unguided rotation-minimising helices up to a quarter turn
with optional twist. Longer helices, curved guides, spline spines and
compound sweeps remain separate coverage work. The P8 threaded-union
refusal remains `intersect.parametric-surface`. [P8 recipe
evidence](../RECIPE.md#work-recorded-for-phase-5)

## Qualification

The 15 ordinary-recipe inputs cover the quarter helix, line guide, twist,
guided twist, reversed motions, a tight-radius helix, and near-contact
and crossing controls. The 90° and 180° guided twists and 30° and 45°
helical twists exercise the edit angles listed in the sweep recipe
evidence. Both tolerances remain 0.01 mm for comparison and 0.005 mm for
the display chord.

Twelve cases pass under 20M work. A circular-profile quarter helix passes
under the 100M allowance, using 49,644,563 work. Its B-rep mesher accounts
for 47,640,841 units. Bringing that display to the 20M target remains
follow-on work. The 13 passed cases cover 345,204 display vertices.
Each verdict, proof, refusal and work receipt matches at 1, 2 and 4 threads
for its allowance. The bound column rounds the certified rational upward
to six decimal places.

| Input | Work at 100M | Composed bound (mm) or refusal |
| --- | ---: | --- |
| Quarter helix | 2,239,094 | 0.008733 |
| Reversed helix | 2,239,504 | 0.008733 |
| Tight-radius helix | 1,247,440 | 0.009992 |
| Circular-profile quarter helix | 49,644,563 | 0.007998 |
| Line guide | 3,964,286 | 0.004064 |
| Twist | 1,824,380 | 0.002859 |
| Guide and twist | 4,006,396 | 0.006240 |
| Reversed guide and twist | 2,135,880 | 0.005563 |
| Line guide with 90° twist | 4,283,702 | 0.009924 |
| Line guide with 180° twist | 8,517,267 | 0.006678 |
| Helix with 30° twist | 2,229,429 | 0.009709 |
| Helix with 45° twist | 4,070,071 | 0.009982 |
| Guide near the spine | 7,147,626 | 0.008244 |
| Guide crosses the spine | 4,252 | `sweep.numerical-definition` |
| Helix profile reaches within 0.01 mm of the axis at the start | 170,695 | `shapes.rational-patch-unproven` |

The two named B-rep refusals occur during fitting. The guide's projected
direction reaches zero. The near-axis helix retains the fitted patch
admission refusal when the projected Jacobian cannot be certified.
