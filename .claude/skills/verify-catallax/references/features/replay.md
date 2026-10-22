# Regression replay

Replay a corpus, resume interrupted work and compare main against a candidate without repeated process startup.

## Sub-features

- `replay.workers`: fresh 1, 2 and 4 worker runs agree.
- `replay.cache`: repeated build/question keys reuse validated artifacts.
- `replay.resume`: complete cases survive interruption.
- `replay.diff`: manifests and strict STEP bytes expose regressions.

## How to get to it (user POV)

Driver commands `run` and `compare`; JSONL and directory inputs are both supported.

## Driving it with control-catallax

Preconditions: built main/candidate binaries, runtime corpus path `INPUT`, external `PROOF` directory.

- Run: `control-catallax.py run "$INPUT" --binary target/release/catallax-verify --workers 1 --cache "$PROOF/cache" --out "$PROOF/one" --json`. Repeat with 2 and 4 workers, distinct output directories and `--fresh`.
- Compare: `control-catallax.py compare "$PROOF/one/run.json" "$PROOF/two/run.json" --json`. Require `equal: true`; repeat for four workers.
- Cache: repeat the first command. Require cache hits and unchanged artifact hashes. Use a different binary to confirm the cache key changes.
- Resume: interrupt a longer run after complete cases appear, then repeat its command. Require complete cases reused and interrupted cases recomputed.
- Candidate: run the same ordered corpus with the candidate binary, then compare the two `run.json` files. A changed STEP byte must fail even if the manifest is unchanged.

## Gotchas

`--fresh` is required for independent load tests. A shared cache proves reuse, not fresh determinism. One command owns each output directory. Case order and duplicate occurrences matter. Nonzero exits and internal errors fail correctness proofs.
