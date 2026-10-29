## 2026-10-25: finite carrier activity and shared literal leaves

Pieces: 1, 5, application benchmarks

Change: rational Boolean construction retains the complete ordered carrier pool.
Existing finite-face exclusion certificates mark which carriers have a
potentially meeting generating face pair. Carriers without one omit the
curve/surface query pass, unless two coincident faces may overlap, which turns
pruning off for that graph. The covered-pair skip still applies between active
carriers. A skipped carrier's crossings with an active carrier come from the
curve/curve query or, when the active carrier covers it, from that carrier's
hits; pairs of two skipped carriers are not queried. When a skipped carrier's
crossing could change the first cut of a kept closed carrier or the order of
isolated points, the graph queries that carrier and restores the complete
pass's order, so section edges and isolated points keep the complete pass's
order. Optional bound witnesses remain in a private model fork, not in the
canonical graph.

The exact literal leaves `-1`, `0` and `1` share immutable storage. Integer
construction avoids a temporary rational for these values. Arithmetic nodes are
not interned or simplified, including expressions with undefined subexpressions.

Why: suffix rebuilds repeatedly solve extended carrier queries that cannot
contribute a finite section, and construct the same small literal leaves during
polynomial and certificate work.

Migration: public function signatures, numeric allowances, geometric tolerances,
STEP options and display controls are unchanged. Separately constructed common
literal leaves may now satisfy `Real::same_node`; node identity is a cache probe,
not a persistent feature name. Optional bounds can add work to an early Boolean.
This is not a claim that every SolidDraft click meets 1000 ms.

The recipe-stage benchmark accepts `--input`, `--repeat`, `--stages` and
`--step-output`. The display benchmark adds `--application` for the real bounded
final-result build, alongside `--incremental` prefix requests. `--edit-from
OLD_NATIVE.json` times a real checkpoint suffix edit and verifies fresh-new and
retained-old OBJ equality outside the timer. Timing, work reporting and scheduling
counterfactuals remain in benchmarks, outside normal part construction.
