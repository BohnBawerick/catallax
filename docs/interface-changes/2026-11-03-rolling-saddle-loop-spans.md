## 2026-11-03: a loop's fitted spans follow its bound

Pieces: shapes
Change: `RollingSaddle::fit` on a loop now halves only the spans over the 100 nm target, instead of doubling them all. `RollingSaddle` gains `span_widths`. `display_grid` now cuts each span by its own width, so its `chord` holds when spans differ. A through ring's spans are all equal, so its fit and its grids do not change.
Why: a loop with one hard stretch paid for it all the way round, and `fit` runs on each rebuild of a fillet step. Without the grid change the chord would have been false on uneven spans.
Migration: code that reads a grid column as an even step in `u` needs updating for a loop. A loop's `span_count` is no longer a power of two times eight.
