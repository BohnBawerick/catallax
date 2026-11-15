# Numerical cone face partitions

`catallax_boolean::numerical::clip_cone_cylinder_face` constructs the retained
regions of one circular cone band outside a cylindrical through-tool. Its inputs
are `&Model`, `FaceId`, a nonempty face key, `&Cylinder`, a common signed axial
interval `[Rational;2]`, `NumericalOptions` and `&mut Decider`. Reuse the original
tool and interval for both cone faces. Coordinates measure millimetres along
the unit tool axis from its original origin.

`None` means an unsupported input face family. The initial family has two
complete circular boundary uses on the positive cone nappe. Supported separated
and point-only cases return `Some`, including their retained material regions.
The tool interval must cover the retained cone source control hull. Blind holes,
apex bands and coincident boundary intervals need further consumers.

`NumericalFacePartition` contains an open `NumericalDefinition`, native cone
trim bindings and all original 3D/UV intersection fits. Every fit retains its
native parameter interval, knots and requested rectangles, including unused
point reports. No native curve is joined, fitted again or reparameterized.

Original rings retain their loop, coedge, edge, both use senses and analytic
circle. `ring_links` identify retained edges, native source angular intervals
and enclosures of the original unwrapped circle angles. These links let the
caller share ring segments with cap and wall assembly.

Each finite event keeps every original native hit and its optional local vertex.
Adjacent searches meeting on a cone angular seam use the ordered roots of the
retained generator against the original cylinder. This avoids treating a small
rounding difference in cylinder UV as a second physical meeting. Root ownership
requires a separated nearest root within the preset physical tolerance. Roots
closer than that tolerance retain distinct IDs. Ambiguous ownership refuses.
Records with the same already-owned vertex are combined without merging other
vertices. Point-only reports never become short edges.

The result includes an immutable `NumericalCylinderSource` and unchanged
tool-side UV templates with canonical quarter IDs. Template binding to final
cylinder wall loop/use addresses belongs to the caller. Source identity is
checked against the finite query's full cylinder source.

Kept arcs and source boundary segments form directed native loops. Finite inward
directions choose the next edge by clockwise order. Sampled loop orientation
and containment assign outer boundaries and holes without an area threshold.
The input face sense is retained. Actual edge/vertex endpoints, native endpoint
coupling, declared endpoints and consecutive physical pcurve joins are checked,
including the closing join. Numerical evidence retains those checks and each
explicit fitted-endpoint move with its before/after position and distance.

The finite producer currently checks its fits at 1e-7 mm. The face result checks
all retained samples against the requested positive tolerance. A tighter preset
can refuse when those fits do not satisfy it. Evidence remains numerical.
There is no exact incidence fact, complete branch-count proof or whole-body
deviation bound. Source and quarter definitions do not authorize STEP export.

The caller owns cap/wall assembly, global material sense, names, edits and final
closed-shell admission through `build_with_native_trims`. `work` records boundary
preparation, producer stages, source/chart admission and topology work under the
caller's allowance and deadline. Standalone calls retain the 20M allowance.

The `numerical_hole_carriers` example accepts `--clip-cone --finite-bands MIN MAX`
to exercise this API on actual E07 prefixes. Reports include the open definition,
source curves, native bindings, ring links, events, wall templates and evidence.

The search reuses the SMLib requested-domain intersector described by the P04
port map. Wire selection follows OCCT's finite inward pcurve directions and
[clockwise edge selection](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/BOPAlgo/BOPAlgo_WireSplitter_1.cxx#L519).
Cone seam ownership follows OCCT's [line substitution into a quadric](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/IntAna/IntAna_IntConicQuad.cxx#L67),
using the existing `RealPoly` quadratic solver and 128-bit root enclosures.
These root definitions identify numerical ownership without changing the
original intersection curves or claiming exact incidence on analytic carriers.
