## 2026-10-28: export-only caps around interior sphere poles

Pieces: 7, STEP consumers

Change: A trimmed spherical face whose parameter pole lies strictly inside it
now joins the partial-chart path of
[the regular periodic trims](2026-10-26-partial-periodic-step-trims.md). An
export-only latitude ring around each such pole splits the face. The pole moves
onto its own cap face, written on a separate `SPHERICAL_SURFACE` whose axis is
the source frame's reference direction, so both chart frames are pole-free. The
holed face keeps its source frame. Kernel topology and geometry are unchanged.

The ring lies on the plane `axis·(x - centre) = lambda·|axis|²` with rational
`lambda`, which keeps the cut free of the axis length's square root. Heights
start at half the radius and move towards the pole. The first height whose closed
cap meets no source trim is used; a certified interior pole then makes the cap lie
strictly inside the face. Faces whose pole is on a trim or uncertified, or whose
rings fail every tried height, keep the whole scene on its legacy encoding, as do
all the other limits of the earlier fragment.

`StepReport.faces` counts the cap faces and `added_edges` includes each ring.
Report fields and manifest formats do not change.

Why: The pinned raw translator needs a degenerate pole edge for a pole inside a
face, and STEP has no such edge. Poles were the most common reason that partial
spherical exports kept their raw-invalid legacy encoding.

Migration: None for callers. Consumers that count faces or surfaces in exported
files see one extra face per interior pole and one extra spherical surface per
split sphere.
