## 2026-10-31: planar draft around one preserved through bore

Pieces: features, field, step

Change: `draft` now answers for a convex planar body with one separated axial
through bore, in both evaluators. The selected planar walls move as before. The
bore wall, its support and both cap rings stay unchanged and keep their names. The
new `catallax_features::draft::bored` module holds the pure proof both evaluators
use: the finite cylinder between two opposed caps stays strictly inside every other
halfspace, on the source envelope and on the drafted one. The B-rep assembles the
drafted cell with the rings. The field drafts the recipe's target and subtracts the
same tool, so it reads no B-rep topology.

`PlanarCell::new` and its `anchored` test now use the supports that define each
vertex instead of proving those exact zeros, and order each face by support
adjacency. This changes work only: vertices, cycles and faces are identical. STEP
decimals now tighten an enclosure before an exact equality proof. The written
digits are unchanged. Field distance no longer proves a foot point on its own plane.

Why: SolidDraft had to suppress a hole before drafting its body, and a rotated
drafted body exhausted its work on zero proofs it did not need.

Migration: callers that matched `draft.unsupported-arrangement` or
`draft.unsupported-source-boundary` for this family now receive an answer. A
selected bore or cap, a bore touching or crossing an outer wall, a blind or
non-axial tool, an oblique neutral plane, several bores, a concave envelope and any
curved outer wall keep those refusals. Recipes, budgets and tolerances are
unchanged. Request work falls for most planar drafts.
