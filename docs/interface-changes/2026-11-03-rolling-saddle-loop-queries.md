## 2026-11-03: rail tables and tested queries on the loop of two cylinders

Pieces: shapes
Change: `RollingSaddle::parameter_table` now answers for a loop instead of refusing. `SaddleLaw` gains `gap_sign`, `narrow_axis`, `side_axis` and `rail_zero`, which gives the direction at a table's angle zero and the way the angle turns, for either kind of ring. `line_hits`, `on_boundary`, `normal_at`, `side_at`, `display_grid` and `range` are now tested on a loop. No recipe reaches a loop yet.
Why: this is the second piece of the slice for curved edges, periodic seams and singular cylinder junctions. A solid on the loop needs these answers before it can classify, display or write.
Migration: code that builds a `SaddleLaw` literal needs `gap_sign`. A caller that relied on `parameter_table` refusing a loop no longer sees that refusal.
