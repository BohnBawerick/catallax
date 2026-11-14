## 2026-11-12: a journal of what work adds to a model

Pieces: 2 (shapes), 3 (topology), 5 (boolean)
Change: `GeomStore` and `Model` gain `journal_begin`, `journal_end` and
`replay`. A journal kept on a copy records every surface, curve, point and
fact, and every vertex, edge, coedge, loop, face with its origins, edge
construction name, blend base, shell, lump and body added to that copy.
`replay` adds the same to the store or model the copy was taken from, in
order, through the same adding methods, and returns a `GeomRemap` or `Remap`
that says where the copy's new handles landed. `Handle::moved` shifts a handle,
and the `Relocate` trait moves the handles inside a result, on both sides of a
`Result`: an error type must implement it too, as `Refusal` does by naming no
handle. `catallax-shapes`
exports `GeomJournal` and `GeomRemap`; `catallax-topology` exports `Journal`,
`Remap` and `Relocate`. In the Boolean, the cell of a face in imprint and the
classification of a piece now run on a copy and are replayed at once. No
answer, refusal, work count or file byte changed.
Why: decision L1 of [PLAN.md](../archive/PLAN-2026-11.md#124-named-for-later). A loop item
that adds to the model must have this form before another thread can run it.
[RUN-AHEAD.md](../RUN-AHEAD.md) section 3 describes it.
Migration: none for callers. A new method that changes a `GeomStore` or a
`Model` must write itself into the journal and be replayed by `replay`. If no
loop item can reach it, it asserts instead that no journal is being kept, as
`Model::map_face_origins` and the writers of the chart expansion tables do.
