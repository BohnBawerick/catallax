# Target parts

These are `catallax-recipe/3` inputs for the roadmap's P0–P6 scoreboard. Dimensions are millimetres, exact decimal or rational strings. P0 is eight separate recipes so an early refusal cannot hide another hard-case family. P1–P6 are complete intended recipes, not reduced substitutes for operations the kernel cannot build yet.

| Input | Fixed dimensions and construction |
| --- | --- |
| `P0/crossing-pipes.json` | Two perpendicular pipes centred at the origin. Outer radius 5, bore radius 3, length 20. Through-bore tools extend 1 past each end. Hollow each pipe, then unite. |
| `P0/line-tangent.json` | Radius-5, height-10 z cylinders whose axes are 10 apart. Union at line contact. |
| `P0/equal-cylinder-sphere.json` | Radius-5 cylinder from z=-10 to 10 and radius-5 sphere at the origin. Union exercises equal radii and tangent contact. |
| `P0/block-sphere.json` | Block from the origin to 10,10,10; radius-6 sphere at 10,5,5. Union. |
| `P0/seam.json` | Radius-5, height-10 cylinders separated by 4 in x. The second uses quaternion 3,0,0,4, an exact rational rotation around z that changes its seam. |
| `P0/near-miss.json` | Same cylinders, separation 10.0000001 in x, a 0.0000001 gap. Same rotated seam. No tolerance-based contact is assumed. |
| `P0/general-cone.json` | Two frusta, radii 5 to 2 over height 10. Second translation 2,1,1 and quaternion 2,1,1,0. Union. |
| `P0/general-torus.json` | Two tori, major radius 8 and tube radius 2. Second translation 1,2,1 and quaternion 2,1,1,0. Union. |
| `P1-bracket.json` | Plate 60×40×5 extruded from a rectangle; flush flange 60×5×25 above z=5. Boss radius 8, height 10, at 30,25,5; through bore radius 4. Four radius-2 bolts at x=10/50, y=12/32. Flange slot spans x=20..40, z=15..21, with radius-3 ends. Eight overlapping drill cuts plus a straight middle cut preserve the rounded slot and exercise the subtraction chain. |
| `P1-drafted-bracket.json` | Complete P1 construction with the flange's `right` wall drafted inward 45 degrees about z=5 before the plate/flange union. The right wall runs from x=60 at z=5 to x=35 at z=30. The boss, bore, bolt holes and rounded slot are unchanged. |
| `P2-collar.json` | Rectangle spanning radii 6..12 and axial positions 0..12, revolved 360 degrees about y. Both outer circular edges have 1×1 chamfers. |
| `P3-filleted-boss.json` | Block 40×30×10, boss radius 8 and height 12 at 20,15,10. Concave junction fillet radius 2, then the four convex top block edges radius 1. |
| `P4-shelled-tray.json` | Centred 40×30×12 block, inward wall thickness 2 and top face open. Intersect with a centred radius-23 cylinder extending from z=-1 to 13 to round its four corners. |
| `P5-flanged-housing.json` | Stepped profile revolved about y. Inner radius 8, body outer radius 12 and height 30; flange radius 22 and thickness 5. A diameter-4 through cutter with diameter-8, 90-degree countersink is circular-patterned four times at 90 degrees on radius 17. Subtract the patterned cutters, then fillet the flange/body junction with radius 2. |
| `P6-lattice-bracket.json` | L bracket within 12×8×8, 3 mm legs and a 1 mm enclosing shell. Cubic infill has 2 mm cells and radius 0.3; its region roots 0.25 into the cavity walls. A radius-0.5 through hole crosses the base. The shell is an explicit subtraction of the documented cavity, so this part does not depend on the separate shell feature. |

Stock quaternions are homogeneous and normalized by the existing exact placement code. Recipes keep operation IDs and origin references fixed. Do not rename them to accommodate a kernel failure. Change dimensions only as a reviewed corpus revision, with a new report and an explanation of why the old input was wrong.

`holes/P1-bolt-holes.json` and `holes/P5-countersinks.json` isolate the named-hole feature gates with the same fixed P1/P5 dimensions. They are not completed target parts: P1's slot chain remains in the full bracket, and P5's pattern and junction fillet remain in the full housing. The hole cases use one placement group rather than replacing the pattern operation's acceptance test.

Run the same report on `parts/holes` to measure these feature gates in both kernels. Application timings still require an actual `catallax-app` build; the scoreboard does not infer them from kernel construction. Paired timing evidence records load averages before and after every AB/BA pair.

For P6 print output, embed the complete recipe in a `catallax-print/1` wrapper with an explicit `tolerance`, then pass that wrapper to the verification driver's `run`. For example:

```sh
jq -n --slurpfile r parts/P6-lattice-bracket.json \
  '{format:"catallax-print/1",recipe:$r[0],tolerance:"0.5"}' > /tmp/P6-print.json
python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py run /tmp/P6-print.json \
  --binary target/release/catallax-verify --out /tmp/P6-print-proof --json
```

The tolerance is an explicit print setting in millimetres, not a change to modelling tolerance or the scoreboard's allowance. An answered request writes `print.stl` and the matching `display.obj`. A printable mesh alone does not complete P6's mass-equivalent STEP stand-in gate.

For a scoreboard that also retains direct field print artifacts, add `--print-tolerance 0.5` to `parts`. Every row then includes a separate `print` response. Unsupported print definitions keep their refusal instead of borrowing the B-rep mesh.

## Run the report

Build `catallax-verify`, then run from the repository root:

```sh
python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py parts parts \
  --binary target/release/catallax-verify --out /tmp/catallax-parts \
  --max-work 1000000 --occt-python /tmp/occt/bin/python --json
```

The optional OpenCASCADE checker runs in a scratch Python environment with `cadquery-ocp`, not in the kernel. Without a binding, each relevant check reports `step.binding-unavailable`. The [exchange contract](../docs/VERIFY-INTERFACE.md#target-part-scoreboard) defines the report and its coverage limits. The workflow runs on `main` and retains the bundle as `target-parts`.

Every recipe runs through B-rep construction, independent field compilation and sampled comparison, even when one refuses. Elementary fillet and chamfer construction and the convex planar shell used by P4 are supported; general contacts and the remaining shell families retain typed unfinished-coverage refusals. P3's display and comparison remain separate mesh-dependent gates. The report does not skip operations or count coverage gaps as impossible geometry. The fixed checkpoint allowance is deliberately below the exchange's default; a work refusal is coverage at that allowance, not a wall-clock timeout or an invalid part.

A passing STEP check means OpenCASCADE opened the file, found valid topology and matched the face count. Imported volume is recorded, but kernel volume comparison is explicitly unmeasured. Application interaction and whole-result tolerance proofs remain separate roadmap gates. None of these recipes counts as a completed target part solely because this report ran.

`edit-checks.json` keeps reproducible parameter edits beside the unchanged target recipes. P0/seam raises cylinder `a` from height 10 to 12, P1 changes the boss radius from 8 to 9, P2 the chamfer setback from 1 to 0.75, and P4 the shell thickness from 2 to 1.5. Each `edit_check` rebuilds both complete recipes, resolves every original face origin and physical-edge construction name against the new incidence, and checks the independent field name registries. Splits and merged origin sets resolve by containment. The seam edit separates previously merged caps and checks that each cylinder rim keeps its own construction pair. Registry lookup is not a field-boundary certificate. Missing names fail; missing features report `waiting` with their actual build refusals. No prefix or replacement geometry counts as the complete target's passing check.

The report checkpoints all configured checks in its top-level `edit_checks` map before starting display or sampled comparison. An unrelated slow display cannot hide the naming results. Each finished geometry row also contains its `edit_check`. Report `status` distinguishes running, completed and interrupted coverage; an interrupted report is not a passing full scoreboard.
