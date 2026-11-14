## 2026-11-11: loops of whole items on every core

Pieces: 1 (numbers), 5 (boolean)
Change: `catallax-numbers` exports `ValuesAhead` and the `values_ahead!`
macro. A loop whose items are values of inputs every thread shares hands its
later items to other threads whole; the request's thread adopts each value
with its charges and statistics in order. The Boolean uses it by default for
probe clearance against boundary edges and for comparisons of a new point
against known points. No signature, outcome, refusal code or work count
changed.
Why: the next step of decision L1 of [PLAN.md](../archive/PLAN-2026-11.md#124-named-for-later).
[RUN-AHEAD.md](../RUN-AHEAD.md) section 5 has the measured times.
Migration: none. An item handed over this way must read nothing a loop
iteration changes and must keep nothing: taking the model by `&Model` is what
shows it. An item that meets a charged cache is not adopted; the request's
thread computes it.
