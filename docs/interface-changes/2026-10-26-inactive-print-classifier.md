## 2026-10-26: inactive definitions in print and mass preparation

Pieces: catallax-field, print and field-mass consumers.
Change: classifier preparation uses the retained plan's existing live-leaf record.
An unsupported historical source leaf no longer causes `field.print-coverage`
when no live instruction references it. Original leaf indices are preserved.
Why: preparation visited definitions absent from the result's evaluation tape.
Migration: no signatures, formats, accuracy, work allowances or output caps
change. Active unsupported definitions still refuse. Consumers must check the
complete print/mass outcome; removing this preparation refusal does not remove
independent refinement limits. See the [behavioral ledger](../INACTIVE-PRINT-CLASSIFIER.md).
