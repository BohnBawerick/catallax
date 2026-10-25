## 2026-10-18: ND3 an empty answer is answered, with no STEP file

Pieces: 7, judge interface, cli
Change: an empty answer has status `answered`, `result.kind` `empty` with zero counts, volume and area exactly 0, centroid `null`, `step` `null` and no `answer.step`. A stale STEP file is removed.
Why: an empty answer can be reached through not-exact decisions (a near miss treated as touching), so it needs the same `exact` flag and marks as any answer, which a separate status would hide.
Migration: the `catallax-cli` creator skips `write_step` for an empty result. Piece 7's `write_step` still writes a valid file for an empty body when called.

