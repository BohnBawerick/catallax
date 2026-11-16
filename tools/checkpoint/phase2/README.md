# Phase-2 scoreboard targets

`checkpoint.py scoreboard` runs these three targets beside P0 to P5 and the census parts, with the same definition of done. Each `*.json` here is a `catallax-recipe/5` file in millimetres. `cases/` pins each recipe's SHA-256 and its one committed scalar edit. CadQuery and SMLib build the same targets from `head-to-head/targets.py` (branch `phase2/`).

They live here, not in `parts/`, so that adding them changes neither the `parts` CI scoreboard nor the published five-tool comparison. The dimensions are this tooling's proposal, written because `docs/TARGET-PARTS.md` left P8 and P9 to be written "when their turn comes" and phase 2 is that turn. The owner can revise them. A changed recipe needs a new hash in its case file.

| Target | Construction | Edit |
| --- | --- | --- |
| P8 helical spring and threaded boss | Plate 50 x 50 x 4 with a radius-8, height-20 boss. A radius-0.6 circle on the boss wall is swept along an 8-turn helix (pitch 1.5) and united: a round thread. A radius-1 circle is swept along a 5-turn helix of radius 17 (pitch 3) that starts buried in the plate, and united: the spring. Both sweeps are rotation-minimising. | Spring rise 15 to 17.5 (pitch 3 to 3.5) |
| P9 lofted duct | A 40 x 20 rounded rectangle (corner radius 5) at z=0 lofted to a radius-12 circle at z=40. The wall is that loft minus an inner loft: a 36 x 16 rounded rectangle (corner radius 3) at z=-1 to a radius-10 circle at z=41. | Inner outlet radius 10 to 10.5, so the round end's wall thins from 2 to 1.5 |
| Curved housing | Smooth loft through circles of radius 10, 14 and 10 at z=0, 15 and 30. Subtract a prism whose floor is a cubic Bezier arch (dipping to z=23.5 at the axis), extruded along y. Fillet the whole cut edge between the loft side and the arch, radius 2. That edge crosses all four quarter patches of the loft side. | Middle section radius 14 to 15, which rebuilds the cut and the fillet |

Why these shapes:

- **P8 profiles are circles in the axial plane.** A plane exactly normal to a helix tangent has no rational quaternion. At these pitches the axial plane is within about 4 degrees of that normal, so the recipe's rotation-minimising frame and OCCT's corrected Frenet frame give nearly the same tube. A triangular thread profile would turn about 23 degrees per turn under a rotation-minimising frame, and the recipe has no screw-motion orientation.
- **P9's wall is a subtraction, not a shell.** A constant-thickness shell of this loft, open at both ends, fails in CadQuery 2.6.1 (`BRep_API: command not done`) and in smlibpy 0.4.0 (`shell failed`, or an exit without a result) at every wall tried, including when the circle is split into eight arcs to match the rectangle. The subtraction builds in both.
- **The curved housing is the chain of roadmap 1.1 and 2.1**: loft, spline cut, fillet across patches, one more edit. Its case asks `check.py` to also open each STEP file with healing off: the repository's raw `StepToTopoDS` builder at 1e-7 mm from `tools/step-regressions/check_numerical_step.py`, with no ShapeFix, sewing or SameParameter. The scoreboard reports that result beside the verdict and does not gate on it. The CadQuery and SMLib exports of P8 and of the curved housing all fail that read while passing the default one; of P9, only CadQuery's passes it. The curved-housing edges carry about 4e-6 mm, which the builder pins to 1e-7 mm. Gating on it would mark both reference kernels not done because of the reader.

No target has an analytic volume or membership formula, so the STEP check decides each one, as for the census parts.

## Census hole position

`cutter.py` keeps the census drill through the centre of the part's two larger bounding-box extents when that removes material. Otherwise it moves the hole to the first of a fixed list of offsets, in eighths of those extents, nearest first. A candidate passes when OCCT's common of the source and the cutter exceeds a thousandth of the source volume at the drilled radius and grows by as much again at the edited radius. The case records every candidate tried. On the ten retained parts this keeps the centre for seven and moves 00000131, 00000632 and 00000633, whose centred hole sat inside an existing bore.
