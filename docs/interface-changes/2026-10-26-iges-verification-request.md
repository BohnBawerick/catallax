## 2026-10-26: explicit IGES verification request

Pieces: judge interface, catallax-cli, catallax-iges consumer

Change: `catallax-verify` accepts source `iges` with problem `catallax-iges/1` (`recipe`, optional `tolerance`, `scale_digits`, `units`). It builds B-rep only, writes `answer.igs` with the existing writer and publishes the `IgesReport` bounds, units, topology counts, byte count and SHA-256 as manifest field `iges`. Writer failures return `iges.*` refusal codes, work exhaustion returns `numbers.operation-budget`, and an I/O failure is a process error rather than a refusal. `about` lists the source and the `iges` format. `answer.igs` joins the artifacts removed before every request, and `control-catallax.py` accepts the wrapper and compares the file bytes. Previously the source returned `cli.invalid-source`.

Why: the merged IGES writer had no replayable request, so the kernel loop could not exercise or compare its output.

Migration: callers that listed the verification sources add `iges`. Existing recipe, STEP and print requests, manifests and bytes are unchanged; their output directories lose a stale `answer.igs` if one is present.
