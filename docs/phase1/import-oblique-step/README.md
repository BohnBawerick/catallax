# Oblique numerical STEP qualification

The fixed sources `00000022` and `00000134` now complete numerical admission,
a pose rebuild, original and edited approximate STEP writes, disk reread checks,
and raw OCCT validation. The source and repaired-definition hashes match the
[admission controls](../import-planar/README.md). Each body has six faces,
twelve edges and eight vertices. Source-to-result upper bounds remain at most
0.000000000001 mm, under the fixed 0.01 mm local tolerance.

Producer: `e89c4c9a012ee78de48e19a4fcbd04d13e3f4788`.
[real-controls.json](real-controls.json) pins the native binary, source,
repaired-definition and actual STEP hashes, charged stages and reference results.

| Public model ID | Complete source/pose/STEP work | Allowance | Raw original and edited STEP |
| --- | ---: | ---: | --- |
| 00000022 | 706,121 | 20,000,000 | Valid |
| 00000134 | 1,111,948 | 20,000,000 | Valid |

The pose translates by 3, 5 and 7 mm. N04's `rebuild` retains the original checked
state; the edited STEP report has one historical stage. Raw reference checks
confirm positive volume, retained topology, volume preservation under the pose,
the translated centre of mass and interior material. Reference centre translation
uses a 1e-8 mm comparison and volume uses 1e-9 times the larger of one and source
volume; these checks do not create a whole-body approximation bound.

The reference uses OCCT Python 8.0.1.1 and the existing raw
`StepToTopoDS_Builder`, with explicit millimetre/radian identity factors, precision
and maximum tolerance 1e-7 mm. It does not use `TransferRoots`, `ShapeFix`, sewing
or `SameParameter`. All input file bytes remain unchanged. The source files and
full geometric reference values remain private under the
[existing licence scope](../../milestones/R8-public-step-census.md#source-and-licence-scope).
This is two real workflows, not a full census or a spline Boolean qualification.

## Reproduction and controls

Build `catallax-import`'s `numerical_import` example in release mode and run it on
each hash-fixed original source, writing to `RESULTS/MODEL_ID`. The example shares
20M work across admission, pose rebuild and both approximate STEP files.
The [writer contract](../../interface-changes/2026-11-15-numerical-step-export.md)
separates numerical source evidence, emitted-file evidence, length rounding and
dimensionless direction-ratio rounding.

Capture the public controls outside the checkout when running these test targets:

```sh
mkdir -p /tmp/catallax-oblique-capture/n03-axis /tmp/catallax-oblique-capture/oblique
CATALLAX_N03_CAPTURE=/tmp/catallax-oblique-capture/n03-axis CATALLAX_OBLIQUE_CAPTURE=/tmp/catallax-oblique-capture/oblique cargo test -p catallax-step --test numerical_step --test numerical_oblique
python3 tools/step-regressions/check_import_oblique.py --real-results RESULTS --controls /tmp/catallax-oblique-capture
```

The reference checker verifies the manifest's source, repaired-definition and
file hashes before checking raw topology, material and pose. The two authored
tubes retain their hole after independent rational rotations, including a
three-component frame whose decimal directions require orthogonalization.
All three original N03 files, including the exact legacy control, remain byte
identical. The exact legacy control retains its 2,400 work units.

Thirteen focused numerical writer controls pass without changing the original
seven N03 tests. Positive UV-roundoff and negative actual plane/endpoint gaps are
separate controls; coarse decimal geometry must pass its own fixed-tolerance
checks. The work-PC Linux format gate, all 270 import tests and all 135 release
STEP tests pass, with eight pre-existing ignored STEP tests. CI-selected
all-target Clippy passes for all eighteen crates selected by the writer and its
new reference checker, including the reverse consumers and numeric dependencies.

Before implementing the extension, inspected OCCT's
[`GeomToStep_MakeVector`](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_9_2/src/GeomToStep/GeomToStep_MakeVector.cxx),
[`StepToGeom::MakeAxis2Placement`](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_9_2/src/StepToGeom/StepToGeom.cxx)
and `gp_Ax2`, plus the [P04 SMLib map](../P04/PORT-CONSUMERS.md).
The writer follows STEP's normalized direction ratios and projected orthogonal
plane frame, then independently checks its selected decimal geometry.
It needs no general projection/tracer port and does not change Blade's curved
trim, spline intersection or Boolean implementations.
