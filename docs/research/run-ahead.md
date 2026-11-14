# Running ahead on other cores in the exact kernel

Outcome: adopted. STEP export is much quicker on many cores. Booleans are not quicker yet.
When: 11 to 14 November 2026.

## What we were trying to achieve

The exact B-rep kernel using every core, with the same answer, the same refusal, the same work count and the same output bytes at any thread count.

## Why it was proposed

Two things on 11 November.

First, a measurement (step R5) of how parallel the exact kernel could be. It counts work units, so it needed no quiet machine. It found that a rebuild of P0 to P5 could run 1.0 to 3.3 times faster in parallel as the code stood, and 1.0 to 8.2 times (median 2.3) with a design called "plan in parallel, append in order". The plan wanted the answer early because every later milestone would otherwise be built on single-thread assumptions.

Second, the owner's ruling the same day: use every core by design, do not gate it on speed studies, keep the rule that the answer is the same at any thread count.

## What was tried

The difficulty is that one request has one work counter and one set of proof caches. Split a loop across threads and the work count starts to depend on the schedule, and the count decides refusals.

Run-ahead keeps one thread, the writer, doing everything in order. Other threads work ahead of it on copies and hand back a result with a tape of the work it cost. The writer adopts the result where it would have done the work itself and replays the tape against its own caches. If no result is ready, it does the work as before.

This went into the STEP export first, then four checks inside the Boolean, then five Boolean loops that add to the model.

## Result

- STEP export of P1: 8.7 s on one thread, 1.1 s on 24.
- P1's Booleans: a median of 2.6 s on one thread and 2.4 s on 24, with all five loops handed over. The workers do deliver results, but the request is no shorter. The note gives two possible causes and says neither is confirmed.
- Work counts, answers and bytes are the same in every run.
- The rest of the Boolean can run ahead behind a setting. It is off by default because it is slower that way today.

One fault showed up at the release: three macOS tests failed on a condition variable used with two locks. It was fixed on 14 November.

## Why it was adopted

It keeps the rule that the machine never changes the answer, it pays at once for export, and the owner had ruled that parallel work is designed in and measured later. Making the Boolean itself quicker is left to the speed phase (milestone 5.1).

## Sources

- [RUN-AHEAD.md](../RUN-AHEAD.md): https://github.com/BohnBawerick/catallax/pull/414 and https://github.com/BohnBawerick/catallax/pull/524
- [EXACT-PARALLELISM.md](../EXACT-PARALLELISM.md): https://github.com/BohnBawerick/catallax/pull/395
- The macOS fix: https://github.com/BohnBawerick/catallax/pull/595
- [archived plan](../archive/PLAN-2026-11.md) section 12.6, row R5
