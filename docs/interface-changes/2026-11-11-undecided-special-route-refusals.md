## 2026-11-11: refusals after an undecided identity test keep its decision

Pieces: 1, 3, 6, 7, 8, 9
Change: `Refusal::keeping_any(Option<NotExact>)` keeps a decision when there is one. `RollingMitreBand::checked` returns the band, or the first exact identity test that stopped undecided while recognising it; `of_face` is unchanged and wraps it. The refusals reached after such a test now keep its decision: the measured fillet ring check, the STEP sphere chart and frame checks, the conical and circular trim matches in features, and the rolling mitre callers in topology, measure, mesh and STEP. Codes, details and control flow are unchanged.
Why: these sites read an undecided test as "not equal", then refused with a piece code, so the refusal lost the numbers code and bound. This is the second evidence group of M12; the first is [2026-11-11-refusal-keeps-undecided-decision.md](2026-11-11-refusal-keeps-undecided-decision.md).
Migration: a caller that refuses when `of_face` returns `None` uses `checked` and `keeping_any`. A search that reads an undecided candidate as a miss keeps the first such decision for the refusal it reaches when nothing matches.
