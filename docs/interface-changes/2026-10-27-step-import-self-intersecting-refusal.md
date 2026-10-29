## 2026-10-27: STEP import refusal `self_intersecting`

Pieces: `catallax-import`, application adapter (SolidDraft)
Change: Adds the refusal `{"code": "self_intersecting", "solid": "#n", "faces": ["#a", "#b"]}` to `ImportRefusal`, and `check_self_intersection(&StepFile)`. `read_step`-based entries (`import_step_bodies`, `build_solid`, `apply`, `rebuild`) return it when two planar faces of one solid intersect, in the source or in a repaired result. It is a geometry claim, like `not_solid`. Rules and limits are in [the contract](../pieces/08-step-import-repair.md#self-intersection). Follows [propose and apply](2026-10-27-step-import-propose-apply.md).
Why: a closed, non-self-intersecting result is required of every import and repair.
Migration: An application that matches `ImportRefusal` exhaustively adds the new case and shows the two face names.
