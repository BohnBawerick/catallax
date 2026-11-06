# Algebraic circle-plane chords

Curve/surface intersection now solves a non-rational circle/plane chord
without eliminating its three world coordinates. Rational definitions keep
their original route. Parallel circle and plane normals keep the existing
section handling. The resulting points still pass the existing branch,
contact, ordering and witness checks.

For circle centre `c`, normal `a`, squared radius `r2`, plane origin `p` and
normal `b`, define `v = a cross b`, `L = v dot v`, `h = (p-c) dot b` and
`D = r2*L - h*h*(a dot a)`. With `L > 0`, the chord midpoint is
`c + (v cross a)*h/L`. Positive `D` gives offsets `v*sqrt(D)/L` on either
side. Zero `D` gives one tangent point. Negative `D` proves separation.
No unit-normal assumption is needed.

The unchanged `acute-5deg` and `acute-halfdeg` recipes now export at 1M work.
Their raw STEP imports have one valid six-face solid each and match independent
OCCT construction volumes. The corpus adds ten material witnesses per case
and checks a separate wedge-section volume formula. The names still denote
approximate angles; dimensions and requested radii are unchanged.

No public API, recipe field, fitted carrier, radius/overflow change, caller
accuracy change, UI change or SolidDraft pin update. Defaults remain 17 digits
and 0.0000001 mm file tolerance. The oblique-volume disagreement and the
remaining construction and consumer gaps are not closed by this work.
