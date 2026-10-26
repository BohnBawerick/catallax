## 2026-10-26: bind complete baseline operation histories

Pieces: verification driver and comparator.
Change: extend [baseline completeness validation](2026-10-26-baseline-archive-completeness.md)
to compare ordered kernel/ID/operation/status histories, including entry counts,
against the independent counterpart retained by the fresh six-run archive.
Why: executed public regressions showed that field checks alone could accept
coherently empty or truncated/renumbered baseline operation arrays.
Migration: identical executables require full non-time log/work equality across
both triplets. Different-build work and decision counters may differ. A changed execution history
requires separate audit before reuse. Legitimate comparison preflights that
execute no operations remain valid with matching independently fresh empty
histories. No producer schema, allowance or SINGLEPAIR behavior changes.
