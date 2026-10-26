# Catallax verification map

Build the release binaries and run the skill's doctor first. Use a separate external proof directory per run. Inputs must be public examples or private runtime paths; never copy private corpora into this repository.

| Feature | Entry points |
| --- | --- |
| [Deterministic answers](answers.md) | Questions, recipes, certified spline bodies and consumers, whole-face definition checks and comparison sources |
| [Regression replay](replay.md) | Persistent workers, cache, resume and strict comparison |
| Qualified printing | `catallax-print/1` wrappers through `run`, with explicit millimetre accuracy and optional `route`. Ordinary parts use B-rep; native-field definitions use field. Inspect `printing.route`, the composed file bound, `print.stl` metadata, matching `display.obj` and actual file hashes. See [the output contract](../../../../../docs/VERIFY-INTERFACE.md#field-printing). |
| Field-mass STEP stand-in | `catallax-stand-in/1` through `run` and source `recipe`, with explicit volume, centroid and inertia Accuracy. Inspect field enclosures, exact surrogate properties, rounded-property inclusion and actual surface/solid import. See [the exchange contract](../../../../../docs/VERIFY-INTERFACE.md#field-mass-stand-ins). |
| [Interactive latency](latency.md) | Pinned matched pairs, old replay and performance report |
| [Operation logs and exported questions](answers.md) | `operation_log`, `operations.json`, portable `question.json`, `replay` |
| [Target-part scoreboard and parameter-edit names](../../../../../parts/README.md) | `parts parts --out DIR`, independent builds, saved P0/seam and P1/P2/P4 edits and OpenCASCADE imports |
| [Released-tool target comparison](../../scripts/head-to-head/README.md) | `head-to-head-reference`, `head-to-head parts`, phase costs, authored P0–P5 edits, independent properties and STEP |

Capture commands, responses and artifact side effects. Report an unavailable entry point as skipped, with its missing precondition. Do not claim another entry point covers it. Cleanup must preserve evidence.
