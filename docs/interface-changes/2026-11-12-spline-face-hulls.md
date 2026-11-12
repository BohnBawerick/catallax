## 2026-11-12: spline face pairs kept apart by their control hulls

Pieces: import
Change: the self-intersection check proves two spline faces of one solid apart before it tries their surface intersection. Faces with no common vertex pass when a plane strictly separates their controls. Faces that share one edge, and no vertex off it, pass when a plane through the edge's two vertices has each face's other controls strictly on its own side, with any controls on the line inside the edge's segment. Only faces whose weights are all positive and whose nets have at most 9 controls are tried. Any other pair is decided as before.
Why: M11 needs a fillet between two spline faces of an imported part, and a part with adjacent spline faces did not import ([prerequisite](2026-11-12-spline-fillet-step.md)).
Migration: none. A solid whose spline faces pass this check used to refuse with `import.self-intersection-not-verified` and now imports. Every other outcome is unchanged.
