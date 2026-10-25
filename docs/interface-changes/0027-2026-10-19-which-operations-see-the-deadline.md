## 2026-10-19: which operations see the deadline

Pieces: all
Change: `ARCHITECTURE.md` now states that only the `Decider` calls that do work (`sign`, `compare`, `enclose`, `real_roots`, `certify_root`) see the deadline, so an operation that reaches its answer without one answers after the deadline. Tests after the deadline accept the `numbers.deadline` refusal or the default-budget answer, unless a description fixes which.
Why: the `Decider` offers no way to read the deadline, so the earlier wording of pieces 4 to 7 promised refusals that a correct build could not give.
Migration: none in code. Test writers follow the rule in their deadline cases.

