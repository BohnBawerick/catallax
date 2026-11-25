Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Tolerance-bearing import topology verification

The [import contract](../../../pieces/08-step-import-repair.md#tolerance-bearing-source-topology) owns the API and persistence semantics.
This capture qualifies storage and tolerance updates before the later healing and body-tier stages.
It adds no new solid admission or census credit.

On 16 Nov 2026 UTC, the Linux timing system producer `69b9d5fd921e5f3268cf30861f0b49fe4ca49dd6` passed Linux `ci/run fmt`, ten `tolerant_topology` tests, twenty existing numerical import tests and six `numerical_cylinder_incidence` tests.
The topology tests exercise shared entities, inconsistent/open declarations, local growth, endpoint propagation, cap boundaries, atomic work exhaustion, closed edges, vertex loops, source units, independent bodies, record tampering and changed-source rejection.
The STEP tests additionally inspect emitted constant V on genuinely oblique parallel and antiparallel circle/cylinder axes and preserve the nonconstant endpoint heights of the tilted skip path.

The [twelve real controls](../../../research/tolerant-import-topology/real-controls.json) reconstruct the supported reader's complete `StepFile` definition exactly and reopen their tolerance records against the same source bytes.
They include all ten retained P06 sources plus `00001413` and `00009338`, with hashes from the frozen 3,000-input selection.
This is a source-topology storage check, not a validity check or a claim that all original STEP information is read.
The actual Cargo artifact was bound to executable SHA-256 `5b08e00cc4a2bbce41342e968ef71d7e32893f45332d091c3c8039e2cc01c139` before running the controls.
The example's header owns its invocation; it writes the tolerance record only after reconstruction and reopen succeed.

Independent OCP 8.0.1.1 checks passed all 256 existing authored circle-incidence STEP files and three material probes each.
All 256 file hashes reproduce the [369-file producer](../step-census-circular-pcurves/README.md).
The [two additional oblique-axis files](../../../research/tolerant-import-topology/oblique-axis-controls.json) each pass builder validity, five-face/nine-edge/six-vertex topology, quarter-cylinder volume within 1e-9 mm³ and three transformed material probes.
The combined check covers 258 actual STEP files and 774 material probes.
Validity includes possible OCCT builder pcurve processing and healing; neither unhealed validity nor a complete source-fidelity bound is established.

The implementation follows the monotone update rule in OCCT V8_0_0 `BRep_TVertex.lxx` and `BRep_TEdge.lxx`, with explicit caller controls, a cap and recorded endpoint propagation.
Their bare C++ constructors start at `RealEpsilon()`; Catallax's stated 1e-7 mm initial value is its own configurable policy, not a claim of constructor-default parity.
`ShapeFix_Root` supplies precision/minimum/maximum controls and `ShapeFix_Shape` dispatches wire/face/shell/solid repair; those later healing operations are not implemented by this storage slice.
The [P04 SMLib map](../../phase1/P04/PORT-CONSUMERS.md) establishes no direct import caller for its three remaining intersection omissions, so this change does not modify those paths.
