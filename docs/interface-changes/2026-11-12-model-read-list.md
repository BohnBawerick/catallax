## 2026-11-12: a list of what work on a copy of a model read

Pieces: 1 (numbers), 2 (shapes), 3 (topology), 5 (boolean)
Change: `FactBook`, `GeomStore` and `Model` gain `watch`, `watched` and
`unchanged_where_read`. Between `watch` and `watched`, a copy lists every
question it is asked whose answer later additions can change: what the fact
book holds and how its classes are tied, how many surfaces, curves and points
there are, the surfaces a closed conic lies on, the coedges of an edge, the
owner of a coedge, the names kept beside faces, edges and bodies, and which
bodies there are. `unchanged_where_read` says whether the store or model the
copy came from still answers each as it did when the copy was taken.
`catallax-numbers` exports `FactRead` and `Watch`, `catallax-shapes` exports
`GeomReads`, and `catallax-topology` exports `Reads`. `ValuesAhead::take_if`
asks the writer whether another thread's value still holds before it takes
the value's charges. In the Boolean, the cut of each face in imprint and the
class of each piece are computed whole by other threads and taken when their
reads hold. No answer, refusal, work count or file byte changed.
Why: decision L1 of [PLAN.md](../PLAN.md#124-named-for-later). A journal made
on a copy is good for the model only while the model still answers what the
work asked as the copy did. [RUN-AHEAD.md](../RUN-AHEAD.md) section 3
describes it.
Migration: none for callers. A new method that reads something of a
`FactBook`, a `GeomStore` or a `Model` that a later addition can change must
write the read into the list and be checked by `unchanged_where_read`. Reading
an entity by its handle needs nothing: what a handle names never changes.
