Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Numerical skew collar hole

`catallax_features::evaluate_numerical_hole` evaluates one terminal through Hole placement and returns a separate `NumericalHoleBuilt`. Callers choose `NumericalOptions` before evaluation. `edit` re-evaluates the authored recipe with explicit options. The original exact evaluator keeps its result type.

The result exposes its recipe, numerical body, source-face origins and actual prefix/assembly work. Original stock faces keep their origin sets. New tool faces retain the compiled Entry/Wall band origins. Stable output keys and native parameter curves belong to the numerical body. The numerical STEP wrapper uses this continuation for supported recipes.

Call `evaluate_with_numerical_hole` when the caller needs ordinary fallback. It returns `HoleEvaluation::Numerical` or `HoleEvaluation::Standard`, retaining the latter's original marks and checks. Unsupported stock resumes the ordinary evaluator from the already evaluated prefix under one shared allowance. Calling the optional entry point and then starting a separate ordinary evaluation would repeat that prefix. Both numerical entry points refuse field stand-ins with `field.stand-in-export` before Hole assembly.

The current continuation accepts an annular cylinder with two outer conical bevels and two planar caps. It requires one constant-radius through tool separated from the bore. Other source families return `None`; unsupported arrangements within this continuation return a refusal. This change does not qualify blind holes, counterbores or a tool crossing the original bore.

## Assembly

The finite Cone/Cylinder and Cylinder/Cylinder producers retain the original compiled Hole axis, reference direction and signed axial interval. The consumer keeps every reported nonconstant wall piece and both native maps. Original ring ownership and requested quarter seams supply endpoint joins. Every join must fit the chosen tolerance.

The cylinder-wall builder consumes directed shared boundaries and adds unmatched quarter seams. It retains all closed regions instead of selecting an area threshold. Every polygon ownership query and every side visited by that query charges work and checks the caller limit. Planar caps use the same shared rings and projected tool sections, with a separate original bore loop. Complete numerical admission checks closure, finite trims, surface composition and joins. Missing exact incidence stays missing.

The approach follows OCCT's connected directed-edge wire assembly and SMLib's requested-domain intersection and curve-on-surface checks recorded in [the P04 port map](../../THIRD-PARTY.md). Shared circular curves can have different rational parameter weights in the cone and cylinder charts. The ring correspondence retains their native parameter and checks the resulting rational map in physical space.

## Finite export parameters

Numerical point samples use an exact rational address with `BSplineCurve::evaluate_box([t,t], ...)`. Converting a decimal endpoint to a binary parameter interval can extend it beyond its finite spline domain. The coordinate result remains an interval; this dispatch adds no exact incidence claim.

The STEP payload extracts each finite 3D spline interval before decimal rounding, as it already does for UV curves. The authoritative body carrier and native edge parameters remain unchanged. The emitted finite basis represents the same interval and passes numerical payload admission under the original tolerance.

This finite embedding also avoids OCCT converting a full closed ring basis to a periodic spline. Its [STEP reader forces periodicity for a closed curve](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/StepToGeom/StepToGeom.cxx#L821). The retained rational ring has unequal end weights; periodic conversion changes its last span. Writing the actual finite edge basis preserves the intended geometry and parameter interval.

STEP remains approximate, with the chosen numerical tolerance and separate scalar/parameter rounding checks. No whole-body error bound or exact result is inferred.

## Persistence and app integration

Parametric persistence saves the original recipe and tolerance, then re-evaluates them on reopen. Cold and live STEP must match. This follows the existing native part's feature-definition persistence model. Numerical B-rep transfer is a separate handover operation with its own cost.

SolidDraft still needs to consume this typed result, retain its numerical category and tolerance, support the authored skew axis, and route preview, names and export accordingly. Kernel recipe replay and raw STEP validation do not establish that app acceptance.
