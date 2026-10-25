## 2026-10-25: one polygonal planar face after a Boolean

Pieces: catallax-features, catallax-field, recipe callers.
Change: the existing single-face thicken fallback also recovers a closed finite
polygonal planar face retained by planar Boolean operations, including a
nonconvex outer wire and holes. Source/offset aliases and every incident rim
origin survive; near selections compare closed trims and keep lost/ambiguous
reference refusals.
Why: the plain-prism and field-leaf restrictions previously rejected retained
planar faces whose normal sweep is still an exact straight prism.
Migration: no recipe fields, dispatch sites, allowances, accuracy or cell caps
change. Keep an explicit single physical face selection (use `near` after a
split). Curved boundaries, multiple selected physical faces and source
provenance not established by this bounded recovery retain typed coverage
refusals. This extends the single-cap contract in
[0165](0165-2026-10-25-sharp-selected-face-thickening.md); it does not qualify
native field display or complete target-part acceptance.
