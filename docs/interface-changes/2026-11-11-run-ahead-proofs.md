## 2026-11-11: the exact kernel's export uses every core; thread count setting

Pieces: 1 (numbers), 2 (shapes), intersect, 5 (boolean), STEP export, verification CLI, judge interface
Change: `catallax-numbers` exports `threads`, `set_threads`,
`set_run_ahead_loops`, `RunAhead`, `AheadSite`, the `run_ahead!` macro,
`Helpers`, `record`, `Charges`, `shared_caches` and `private_cache`. The STEP
pcurve fit hands later spans to helper threads through them, on every
available core by default. The Boolean's loops hand later items to workers
only when `CATALLAX_RUN_AHEAD=loops` or `set_run_ahead_loops(true)` asks,
because that is slower today. `catallax-verify` and `catallax` take
`--threads N` before the command, and every program that links the kernel
reads `CATALLAX_THREADS`. Helper threads work to the writer's deadline, and
`Charges::adopt` takes the writer's decider and adopts nothing past it. `uv::fit` in the STEP crate now needs `Sync` pieces,
so the pcurve placement and surface are held in `Arc`. No operation signature,
outcome, refusal code, budget default or file format changed.
Why: decision L1 of [PLAN.md](../archive/PLAN-2026-11.md#124-named-for-later), adopted on
11 November 2026. [RUN-AHEAD.md](../RUN-AHEAD.md) describes the mechanism and
the measured times.
Migration: none for callers. Answers, refusals, work counts, statistics and
file bytes are those of one thread at any thread count and with either
setting. A cache on a shared object, or one a thread keeps for itself, whose
miss spends decision work must ask `catallax_numbers::shared_caches()` before
it is read or filled, or call `private_cache()` when it is read, and must be
listed in `crates/catallax-boolean/tests/shared_state.rs`. A loop of
independent proofs can be marked with `run_ahead!`; the item body must own
what it uses and run on its own copy of the model.
