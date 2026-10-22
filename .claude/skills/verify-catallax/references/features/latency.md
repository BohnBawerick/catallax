# Interactive latency

Measure latency separately from geometric correctness and retain slow cases for the kernel lanes.

## Sub-features

- `latency.pairs`: pinned, alternating AB/BA measurements.
- `latency.legacy`: original process-per-question elapsed replay.
- `latency.report`: target violations grouped by operation, shapes and outcome.

## How to get to it (user POV)

Driver commands `pair`, `legacy` and `report`.

## Driving it with control-catallax

Preconditions: doctor identified an allowed `CORE`, both stream binaries exist, external `PROOF` directory.

- Pair: `control-catallax.py pair "$INPUT" --baseline "$BASE" --candidate "$CANDIDATE" --cores "$CORE" --out "$PROOF/pairs" --json`. Require one warmup and three measured pairs per case in `pairs.json`.
- Legacy: `control-catallax.py legacy "$INPUT" --binary "$LEGACY" --cores "$CORE" --out "$PROOF/legacy" --json`. Capture total elapsed time and each exit status.
- Report: `control-catallax.py report "$PROOF/one/run.json" --json`. Inspect `outcomes` separately from `slow`. Targets are 100 ms stock, 500 ms boolean, 1,000 ms recipe and 2,000 ms full comparison.

## Gotchas

Use identical inputs and allowances. Record ambient load and affinity. Cached wall times cannot establish a fresh latency result. A faster refusal is not evidence that geometry became faster. Group labels suggest profiling targets; they do not identify measured hot functions.
