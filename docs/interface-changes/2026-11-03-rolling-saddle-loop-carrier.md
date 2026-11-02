## 2026-11-03: fitted ball carrier on the one loop of two overlapping cylinders

Pieces: shapes
Change: `RollingSaddleDefinition` gains `overlapping`, for two nonparallel cylinders whose shadows along the narrow axis cross, and `RollingSaddle::fit` fits that loop as a closed rational carrier of degree (12, 2), with rails exactly on both cylinders, a recomputed two-sided bound and a regularity proof. `SaddleRing` is new and `SaddleLaw` gains `ring`. `RollingSaddle::parameter_table` refuses a loop with `shapes.rolling-saddle`. No recipe reaches a loop yet.
Why: this is the first piece of the slice for curved edges, periodic seams and singular cylinder junctions. Both cylinders' angles turn round on this loop, so the chart of a through ring cannot follow it, and it needs a chart split at those critical points or one with none.
Migration: code that builds a `SaddleLaw` literal needs the new field. Nothing else that compiled before changes.
