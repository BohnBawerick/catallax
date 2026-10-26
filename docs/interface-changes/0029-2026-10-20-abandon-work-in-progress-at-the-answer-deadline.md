## 2026-10-20: abandon work in progress at the answer deadline

Pieces: 1, 4, 5, 6, 7, cli, judge interface
Change: decisions check the deadline during budgeted work, not only on entry. `Decider::check_deadline()` provides a checkpoint without changing decision statistics. `TracePiece::{evaluate, derivative, taylor}` now take `&Decider` and return `Result<_, Refusal>` so series evaluation and subdivision can stop. `Ball::{pi_checked, atan_checked}` share the unbudgeted algorithms but check the deadline during their series. Measurement and STEP fitting propagate these refusals without retaining partial results. The judge limits paragraph now promises a `numbers.deadline` manifest and exit 0 within 1 s of the deadline, including work in progress.
Why: release cylinder-sphere questions exceeded a 2 s deadline by 4.751 s, and some exceeded a 50 s deadline without a manifest at 75 s. Checking only the start of a decision did not protect the 60 s outer limit. This supersedes the 2026-10-19 entry about which operations see the deadline.
Migration: pass the current decider to trace enclosures and propagate failures. Long loops outside numeric decisions can call `check_deadline()`. Use the checked pi and arctangent enclosures inside deadline-bearing operations. File formats, deterministic budgets and the determinism exception for deadline refusals do not change.

