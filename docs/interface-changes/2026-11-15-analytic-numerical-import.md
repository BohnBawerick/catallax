# Finite analytic numerical source admission

`catallax_import::numerical::import_step` keeps successful certified imports on
their existing path. A failed planar repair can now retain shared source 3D
geometry and reconstruct finite parameter trims for planes, cylinders, lines and
circular arcs. Missing repair certificates can take this path; a proven violation
of the source cap or neighbouring-feature limit remains a refusal.

`NumericalImport::source_admission()` returns an optional
`NumericalSourceAdmission` with the actual strict refusal, fixed tolerance and
versioned reconstruction method. The method
`retained_shared_analytic_boundaries/1` preserves source vertices and carriers.
Its source record contains no invented whole-face bound or derivative hash.
Body evidence lists actual structure, orientation, finite-loop and residual
checks; its empty movement history means no 3D vertex was moved. The reconstructed
parameters and source interpretation remain separate from geometric movement.

N01's numerical definition now accepts finite circular arcs and finite
cylindrical charts. A circular arc must span less than one turn. A cylindrical
parameter polygon must remain within one finite, non-overlapping periodic chart;
noncontractible complete bands need an explicit seam construction and refuse.
Plane/circle loop checks use sampled projected walks, as planar native splines
already do. These checks establish neither continuous between-sample separation
nor whole-body self-intersection freedom or a complete approximation bound.

The approximate STEP dispatcher adds an analytic writer. It retains explicit
shared edges/vertices and cylindrical pcurves, reconstructs the actual rounded
STEP frames and rechecks decimal 3D geometry. Separately recorded samples check
the actual decimal pcurve mappings. Pcurve parameter rounding reports lexical
numeric selector error across coefficients with different units; this is no
common-unit geometric bound. Actual 3D distance samples enforce tolerance. It
stays separate from length-scalar and 3D direction-ratio errors.
Files remain marked numerical approximate without a global uncertainty claim.
The existing finite planar writer and exact exporter keep their paths.

No exact model extraction or exact-incidence fact is added. The unchanged 20M
request allowance covers admission and every continuation/export stage. Source
bytes are charged before the initial verification hash, including a mismatch;
the already verified digest is reused by subsequent source records. Caller
recipes, source bytes and previously checked bodies are not changed on failure.
The SolidDraft numerical-group and replacement UI integration remains separate.

See [the analytic qualification](../phase1/import-analytic-numerical/README.md).
