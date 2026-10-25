## 2026-10-18: ND9 running environment, limits and exit codes

Pieces: judge interface, cli
Change: the default `deadline_ms` becomes 50,000, counted from the program's start. The judge's outer limits are 60 s and 4 GiB; the kernel does not measure its memory. The program needs an empty sandbox only. Reruns of one build are byte-identical except `stats` and deadline refusals. A refusal exits 0; every nonzero exit, signal, timeout or missing manifest is a crash.
Why: a 50 s deadline leaves time to write the refusal manifest inside the judge's 60 s limit, which a 60 s deadline would race.
Migration: the `catallax-cli` creator applies the new default and sets the decider's deadline from the program's start.

