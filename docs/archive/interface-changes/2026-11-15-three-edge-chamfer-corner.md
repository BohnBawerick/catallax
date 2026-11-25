Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Three-edge planar chamfer corners

Three convex line bevels at a regular three-plane vertex now meet in a sewn
planar triangle. Each input face supplies the intersection of its two contact
rails. Those three vertices trim the bevel strips and define the corner plane.
The original vertex determines the removed-material side and outward winding.
Each triangle edge is shared with one strip through the existing boundary
replacement code. The existing four-top-strip sewing path is unchanged.

The corner retains a `Role::Corner` origin from the input vertex name.
Each strip retains its ordered setbacks and `Role::Bevel` origin.
The implementation adds no recipe format or application API.
The model and its line-plane contacts remain exact where the existing cheap
plane predicates settle them.

The construction follows the planar case of OCCT's
[`ChFi3d_ChBuilder::PerformThreeCorner`](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder_C3.cxx).
The [P04 SMLib map](../phase1/P04/PORT-CONSUMERS.md) supplies no missing tracer
consumer for these existing plane/rail intersections.
No spline intersection or general curved corner filling is introduced.

The public fixture is
`crates/catallax-features/tests/fixtures/chamfer-three-edge.json`.
`crates/catallax-features/tests/chamfers.rs` checks equal setbacks, ordered
unequal setbacks and a setback edit. All three controls match the OCCT
7.9.3.1.1 reference volumes, have ten planar faces, retain the same face
origins, classify material on both sides, and export STEP.
The existing four-top-strip regression remains the sewing control.

Field corner construction retains its previous refusal.
SolidDraft editor qualification and independent raw STEP checks are recorded
with the delivery evidence; this interface change alone does not close the
milestone's app acceptance gate.
