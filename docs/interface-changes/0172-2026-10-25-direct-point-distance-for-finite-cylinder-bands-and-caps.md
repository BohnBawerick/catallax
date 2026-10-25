## 2026-10-25: direct point distance for finite cylinder bands and caps

Pieces: 3, 6

Change: `catallax_shapes::distance` adds `Finite` and `Branch`, an exact recognizer
and direct point-distance reduction for complete circular cylinder bands and
perpendicular circular caps. `catallax_measure::distance::PreparedFace::point_to_face`
uses it for recognized faces only.

Why: repeated chart subdivision dominated explicit and selector distance queries
on these faces although the nearest point has a closed form.

Migration: none. Existing signatures, options, tolerances and the complete search
for every other trim are unchanged. Recognized faces now return tighter bounds with
`cells == 0`; both bounds remain valid.

