## 2026-10-23: exact elementary fillet construction

Pieces: features, field, comparison

Change: recipe fillets now construct exact cylinder strips on plane-plane line edges and torus strips on perpendicular plane-cylinder rings. Shared edge replacement assembly trims inherited loops, inserts the strip, sews contact and transversal-end curves and retains face origins. `Built::fillets` stores each `FilletDefinition` with source supports, orientation, input edge, radius, material side and exact contact curves. The field independently resolves contact intervals and material sectors, compiles the edge regions and exposes the same radius-independent blend origins. Comparison no longer skips fillet operations as wholly unsupported.

Why: these elementary blends have exact analytic carriers and need no sphere-union approximation or fitted spline face. Construction and export must remain separate from display readiness.

Migration: no recipe schema, command, budget, tolerance, legacy question or STEP-option change. Handle the specific `fillet.*` and `field.fillet-contact` construction refusals for unfinished domains. General contacts, roll-over, corner patches, filtered field references and nonrational sections remain uncovered. P3 application display and whole-part comparison are not acceptance claims of this change; both still require the periodic trimmed-face mesh prerequisite. SolidDraft's positional-selector converter is unchanged.

