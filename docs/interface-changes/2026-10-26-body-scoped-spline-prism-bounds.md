## 2026-10-26: spline-prism bounds respect cross-body entity sharing

Pieces: measurement and application callers
Change: body bounding queries derive spline-prism adjacency from face uses in
the requested body. Coedges owned by another body, or not yet attached to any
face, cannot enlarge that body's complete finite-domain inventory. This repairs
an acceptance gap in [the original bounds outcome](2026-10-26-certified-spline-prism-bounds.md).
Why: topology explicitly permits different bodies to share faces, loops,
coedges and edges. Model-wide first-owner back-links are not body ownership.
Two valid prisms sharing one physical cap ring were incorrectly refused despite
both passing `Model::check` and `SplinePrism::check`; no false enclosure was
returned.
Migration: public signatures and controls are unchanged. Each private query
fixes its immutable model, optional body scope and requested `Accuracy` for the
lifetime of its cached extrema. Face-only queries retain their model-wide
semantics and explicitly refuse ambiguous prism adjacency with the existing
`measure.unsupported-surface` code. Every accepted component still passes the
complete prism checker and the unchanged extremum-width proof. Missing caps
cannot be borrowed from another body. Model bounds remain separate from exact
target bounds and their certificate composition.
