## 2026-10-26: explicit qualified baseline evidence reuse

Pieces: verification driver and comparator.
Change: add opt-in `qualify`, complete baseline build-evidence keys and explicit
directory archives of earlier fresh six-run strict qualification. Candidate
1/2/4 always executes fresh; baseline selection revalidates archived evidence
with the current comparator. Run reports now retain fresh policy and worker
exit codes. Comparison accepts artifact paths relative to each run report.
Why: make bounded baseline reuse reviewable without accepting unknown identity,
cached candidate work, stale equality flags or incomplete proof artifacts.
Migration: existing `run`, `compare`, `pair` and parts callers retain their
interfaces. Old evidence lacks required provenance/exits and cannot seed reuse.
Use [the operational contract](../../.claude/skills/verify-catallax/scripts/baseline-reuse.md)
for key fields, strict preservation limits, audit behavior and timing scope.
Parent-owned real-corpus FULL SIX remains required after green CI for this
driver/comparator change. No allowance, geometry tolerance, changed-export
audit or determinism gate is weakened.
