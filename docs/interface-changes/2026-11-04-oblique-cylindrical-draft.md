# Oblique cylindrical draft

The existing recipe `draft` operation now builds a nonzero oblique-neutral cylindrical cell when the selected complete cylinder is parallel to pull, both complete cap planes are perpendicular to pull, and every section has the proved positive convex-section margin in RECIPE section 4.10. No recipe version, operation parameter or tolerance setting changes.

## SolidDraft library consumers

- The wall is `Surface::BSpline(SplineSupport::ObliqueDraft(_))`, with a retained exact ruled law and a rational degree `(4,1)` carrier. It is not a fitted surface or a cone. The source wall and cap origins remain selectable.
- The model has one side, two caps, two physical spline rings and no artificial seam vertices. `DraftedCylinder::check` checks the complete boundary and its material orientation before geometric consumers use the finite law.
- B-rep and native-field display share the checked geometric definition, not topology or a B-rep fallback. Each publishes a closed grid with the requested chord bound. B-rep preview retains pickable face/edge handles and carries no certificate. A certified angular control still returns `mesh.oblique-draft-angular-coverage`; callers must not advertise an angular certificate from preview sampling.
- Requested B-rep volume is available without computing area, centroid or inertia. Those other requested properties keep their existing typed coverage refusals. Field point classification is exact; its current nonboundary distance enclosure is conservative and must not be treated as a tight proximity certificate. The tentative fast scalar is a radial/cap residual, not true distance.
- STEP splits the one side into two open rational charts. Shared rounded cap rows preserve planar incidence, and a checked radial homotopy proves the output boundary and its separate file-to-model displacement. `StepReport::face_map` maps the two side patches back to the same model face. File topology has four faces, six edges and four vertices. No modelling tolerance is consumed. The writer requires rational direction frames; an irrational placed frame retains `step.spline-proof` rather than rounding a new modelling frame.

Hollow oblique arrangements, general cap trimming, swallowed oblique caps, interacting curved components and unresolved remote neutral-plane retention remain typed refusals. `keep_regions` still selects strict material interiors, including when the neutral plane lies outside the body. An explicit empty list selects the empty result and retains consumed-face provenance.

## Verification exchange

An ordinary recipe producing this exact ruled body returns `catallax-manifest/5`. Its record has the manifest/3 recipe fields, but permits the genuine `bspline` surface/curve tags in `result`. It is not a manifest/4 tolerance certificate and does not invent certificate or exception records. `about.exact_ruled_manifest` advertises the version. Existing analytic recipes and refusals retain manifest/3; fitted certified recipes retain manifest/4. Consumers must inspect the format before interpreting the geometry vocabulary.

The public fixture is `crates/catallax-features/tests/fixtures/draft-oblique.json`. The raw STEP regression suite is `tools/step-regressions/test_draft_oblique.py`; set `CATALLAX_VERIFY_BINARY` for the built verifier and `CATALLAX_DRAFT_PROOF` to retain its recipes, driver responses, actual STEP/OBJ artifacts and independent polar-volume observations outside the repository.
