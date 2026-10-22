# Catallax verification map

Build the release binaries and run the skill's doctor first. Use a separate external proof directory per run. Inputs must be public examples or private runtime paths; never copy private corpora into this repository.

| Feature | Entry points |
| --- | --- |
| [Deterministic answers](answers.md) | Questions, recipes, whole-face definition checks and comparison sources |
| [Regression replay](replay.md) | Persistent workers, cache, resume and strict comparison |
| [Interactive latency](latency.md) | Pinned matched pairs, old replay and performance report |
| [Operation logs and exported questions](answers.md) | `operation_log`, `operations.json`, portable `question.json`, `replay` |
| [Target-part scoreboard](../../../../../parts/README.md) | `parts parts --out DIR`, independent builds and OpenCASCADE imports |

Capture commands, responses and artifact side effects. Report an unavailable entry point as skipped, with its missing precondition. Do not claim another entry point covers it. Cleanup must preserve evidence.
