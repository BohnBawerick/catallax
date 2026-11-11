## 2026-11-11: work-span observer behind a development feature

Pieces: 1 (numbers), 2 (shapes), 3 (topology), intersect, boolean, features, verification CLI
Change: a new off-by-default `work-span` feature on `catallax-numbers`, forwarded
by `catallax-shapes`, `catallax-topology`, `catallax-intersect`,
`catallax-boolean`, `catallax-features` and `catallax-cli`. With it,
`catallax_numbers::work_span` exposes `with_work_span`,
`with_work_span_on_this_thread`, `region`, `steps`,
`parallel`, `search`, `append`, `read`, `write`, `write_on_drop` and `key`,
and its report types. Kernel loops, arenas, the fact book, request proof caches
and Boolean caches call these hooks only under the feature. The `work_span`
example in `catallax-cli` answers one request and reports its work and span
([EXACT-PARALLELISM.md](../EXACT-PARALLELISM.md)).
Why: PLAN R5 measures how parallel the exact B-rep could be before decision L1.
Migration: none. Default builds compile none of it and charge the same work.
Observation never charges work or changes an outcome. Code that adds a parallel
or search loop, a shared cache or a new arena-like store should add the matching
hook under the feature so later measurements keep their dependencies.
