# NURBS toolkit

Milestone 2.2 constructs exact rational surface definitions. It does not create sewn bodies or certify regularity or injectivity. Those remain topology decisions. The owner's approval to adapt SMLib surface-building code widens the 10 November port rule, recorded here on 17 November 2026.

## Public API

`catallax_shapes::spline_tools::RationalSplineCurve` has `insert_knot`, `raise_degree`, `split`, `joined`, `reversed`, `reparameterized` and exact rational `point_at`. Knot insertion and splits retain the original parameters. A join gives each directed piece one unit interval and requires exactly matching endpoints. Degree elevation keeps each original span as a Bezier span. The geometry retains its continuity, but the resulting representation only promises C0 internal knots.

`RationalSplineSurface` has `through_grid`, `coons`, `gordon`, `guided_skin`, `n_sided_fill`, tensor-direction `insert_knot`, `raise_degree`, `split` and exact rational `point_at`. Controls and grid points are u-major, v-minor.

- Grid interpolation takes explicit increasing u/v parameters and degrees. Averaged knots and rational linear solves interpolate every supplied point.
- Coons takes two u-running boundaries at v-min/v-max and two v-running boundaries at u-min/u-max. Domains map affinely to [0,1]. Corners must agree exactly.
- Gordon takes u-running profiles at explicit v parameters and v-running guides at explicit u parameters. Both skins use cubic interpolation, or lower degree for fewer curves. All supplied intersections must agree exactly. The result is the sum of the skins minus the interpolating intersection grid.
- Guided skin applies that Gordon construction to cross sections and one or more guide curves. A single guide is supported. The guide must meet each profile at the stated parameters.
- N-sided fill takes a directed closed boundary and an explicit centre. Each boundary splits at its parameter midpoint. Shared straight spokes and Coons quadrilaterals produce N patches with C0 joins. This does not claim G1 continuity or silently sew the result into a body.

No input curve moves, no knot snaps to a nearby value, and no numerical fit is labelled exact. Rational Gordon networks use Bernstein numerator/denominator products instead of replacing rational curves with polynomial approximations.

## Reference methods

The donor pin is SMLib `c9c5979f04cea8c38213004d4f04df0dcfb40d05`. Adapted units retain Apache-2.0 notices and appear in [the inventory](THIRD-PARTY.md).

| Algorithm | References checked | Method used |
| --- | --- | --- |
| Knot insertion, split, elevation, join | SMLib [NL_CrvTool.cpp](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/source/SMLib/src/NL_CrvTool.cpp), OCCT [Geom_BSplineCurve](https://occt3d.com/dev/doc/refman/html/class_geom___b_spline_curve.html) | SMLib homogeneous insertion and Bezier elevation, exact rational arithmetic. Extraction supplies split endpoints. Join checks exact endpoints and scales whole-piece homogeneous weights. No donor parameter snapping or numerical knot removal. |
| Grid interpolation | SMLib [NL_SrfFit.cpp](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/source/SMLib/src/NL_SrfFit.cpp), OCCT [PointsToBSplineSurface](https://occt3d.com/dev/doc/refman/html/class_geom_a_p_i___points_to_b_spline_surface.html) | SMLib tensor interpolation with explicit parameters, averaged knots and exact solves. No least-squares approximation. |
| Coons, Gordon, guided skin | SMLib [NL_SrfAdv.cpp](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/source/SMLib/src/NL_SrfAdv.cpp), [SmBSplineSurfaceNLib.cpp](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/source/SMLib/src/SmBSplineSurfaceNLib.cpp), OCCT [GeomFill_Coons](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKGeomAlgo/GeomFill/GeomFill_Coons.cxx) | SMLib skin-plus-skin-minus-grid construction. Coons is the two-boundary case. Rational arithmetic replaces polynomial approximation, corner averaging and curve movement. OCCT code is not copied or adapted. |
| N-sided fill | SMLib `CreateCornerBlend` in the same surface unit, OCCT [BRepFill_Filling](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKBool/BRepFill/BRepFill_Filling.hxx) | Explicit C0 Coons decomposition with shared spokes. Neither OCCT plate deformation nor SMLib's G1 cross-tangent solver is copied. |

## STEP delivery

`catallax_step::write_spline_surfaces_step` writes an AP242 geometric surface set, not a sewn shell or solid. Each nonempty knot rectangle becomes one Bezier surface with exact [0,1] knots. This prevents decimal knot rounding from introducing an unmeasured parameter displacement.

The caller states a positive tolerance before export. Decimal rounding is approximate unless every emitted control and weight is unchanged. Both the file description and `SplineSurfaceStepReport.exact` record that distinction. The report records the requested tolerance and a whole-surface distance bound. Excessive rounding or a nonpositive rounded weight refuses before writing bytes.

For each patch, let `e` bound the L1 control-point rounding error, `d` the weight rounding error, `r` the greatest L1 distance of a source control from the first control, and `w` the smallest emitted weight. Positive Bernstein bases give the bound `e + 2*r*d/w`. It bounds Euclidean distance as well. The largest patch bound is the file bound. Exact construction does not imply exact decimal STEP.

## Verification

On the work PC through `cx-offload`:

```sh
cargo run -p catallax-step --example nurbs_toolkit -- OUTPUT_DIRECTORY
```

The runnable example checks nonuniform point-grid interpolation, a 3-by-3 Gordon network, single-guide profile preservation, rational Coons boundaries with different degrees and knot vectors, nonplanar five-sided fill joins, weighted curve insertion/elevation/split/join, and surface insertion/elevation/split. It writes six STEP artifacts. The recorded run used 16,461 work units.

The six files opened with `cadquery-ocp 8.0.1.1.0` through direct `StepToTopoDS_Builder.Init` on their geometric sets. No `TransferRoots` or `ShapeFix` ran. `read.stdsameparameter.mode` and `read.maxprecision.mode` were zero. All 16 emitted faces passed `BRepCheck_Analyzer`. Independent OCCT evaluation of the grid, Coons, Gordon and guided-skin laws had maximum sampled L1 error below `3.2e-15 mm`. This sample check is separate from the writer's whole-surface rounding bound. The two approximate exports had bounds below `1.4e-15 mm`, against a tolerance fixed at `1e-6 mm`.

Surface-body operations, G1 N-sided joins and app commands are not established by this surface-definition and exchange proof. The roadmap's next surface-body milestone owns sewing, trimming, extension and thickening.
