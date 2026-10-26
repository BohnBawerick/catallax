## 2026-10-26: bounds for checked finite spline prisms

Pieces: measurement and application callers
Change: the existing `bounding_box` and `face_bounding_box` functions answer for
the represented inward elliptic-offset spline-prism family admitted by
`SplinePrism::check`. Both complete planar caps and the finite ruled side use
their actual shared closed spline rims. Each minimum and maximum must satisfy
the requested `Accuracy` before an answer is published. Signatures and formats
are unchanged; the former `measure.unbounded-face` refusal on this family is
repaired.
Why: an endpoint-free closed spline rim is finite, and a checked active carrier
domain is distinct from an arbitrary face trim or the original offset law.
Migration: existing callers use the same functions and immutable model/face
handles. These boxes enclose the represented B-spline geometry without display
padding or fitting-error padding. For the retained exact target, explicitly
compose the existing complete correspondence certificate: a uniform bound E
allows expanding every coordinate endpoint by E, with a correspondingly wider
extremum enclosure. That composed width is not the original requested model
accuracy. Arbitrary unproved spline trims and geometry declined by the complete
prism checker still refuse. General spline extrema, intersections and new
carriers are not introduced. See the [bounds ledger](../SPLINE-BOUNDS-LEDGER.md)
for the finite-domain proof, independent oracle and remaining limits.
