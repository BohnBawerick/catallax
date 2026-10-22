# Deterministic verification exchange

`catallax-verify` is an opt-in command. `exact-brep` and `catallax` retain every behavior in [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md), including their wall-clock deadlines. This contract does not replace that exchange.

## Commands and framing

Build with `cargo build --release -p catallax-cli --bins`. `catallax-verify about` prints protocol versions, sources, kernel identity and the default work allowance. `catallax-verify stream` reads one UTF-8 JSON object per line on stdin and flushes one response per line on stdout. EOF shuts down the process. Each line is a complete independent request. Diagnostics go to stderr. A worker processes requests in order; independent workers may process different requests concurrently.

```json
{"format":"catallax-request/1","source":"question","problem":{"format":"exact-brep-question/2","id":"box","operands":{"a":{"shape":"block","size":["1","2","3"]}},"steps":[],"result":"a"},"out":"/tmp/box-proof","max_work":100000000}
```

All fields except `max_work` are required. Unknown envelope fields are invalid. `out` names an owned artifact directory; the command replaces `manifest.json` and `answer.step` there. Never assign the same output directory to concurrent requests. `max_work` is an unsigned 64-bit integer; zero is valid and refuses on the first work checkpoint. The default is 100,000,000.

Sources:

| Source | Problem format | Result |
| --- | --- | --- |
| `question` | `exact-brep-question/2` | Full geometric manifest and optional STEP. Parse all legacy fields, but ignore `limits.deadline_ms` during evaluation. |
| `recipe` | `catallax-recipe/3` | Build through `catallax-features`, summarize topology and export STEP with default export options. Recipe properties are not measured by this adapter. |
| `compare` | `catallax-recipe/3` | Run `catallax-compare` with its default options. Report sampled `agree`, `disagree` or `unresolved`, findings and unsupported operations. No STEP output. |

A new source implements one evaluation adapter returning a manifest. External corpus adapters export Python `load(Path)` and yield `{"source":...,"problem":...}` records. These adapters are trusted executable code selected explicitly with `--adapter`. Data files cannot select executable adapters.

## Response and geometric manifest

A `catallax-response/1` object has `manifest`, `work` and `wall_ms`. `work` contains `used`, `limit` and `exhausted`. `wall_ms` is measured elapsed time and never changes an outcome. Work counts are deterministic for the same kernel build, input and allowance. They may change when the implementation changes.

The geometric manifest is `catallax-manifest/3`. For question sources it retains the fields and meanings of manifest/2. Recipe sources return `status`, `exact`, `result`, `marks`, `step`, `refusal` and `kernel`; omitted properties mean unmeasured. A comparison returns `comparison` instead of geometry. `comparison` contains `verdict`, sample counts, ordered findings and `not_comparable`. Findings contain `code`, `detail`, `pass` and diagnostic exact-coordinate `witness`, `brep_faces`, `field_faces` and `gap` strings. These are sampled display checks, never global geometric certificates; `exact` is false for this source. Free-text details and witness formatting are diagnostic.

The surface vocabulary is `plane`, `cylinder`, `elliptic-cylinder`, `cone`, `sphere`, `torus`. The curve vocabulary is `line`, `circle`, `ellipse`, `parabola`, `hyperbola`, `intersection`. Elliptic cylinders use their own tag when their builders land. Parabolas and hyperbolas retain their own tags and finite trimming topology. No old tag is remapped. New analytic kinds require another version change; unknown library kinds must not be published under a convenient existing tag.

`numbers.operation-budget` is a whole-request refusal. It clears partial geometry, marks and STEP output. No deadline refusal originates in this command. Structural input errors have status `invalid-question`; a panic has status `internal-error`. A worker I/O error exits nonzero and is an infrastructure failure. External termination or machine resource exhaustion leaves an incomplete request, never a geometric refusal. A verifier must fail correctness validation for an internal error, missing response or nonzero exit, even if another build fails the same way.

## Work policy

Each `Decider` and `Work` cancellation checkpoint costs one unit. Accepted arithmetic reservations also cost their existing deterministic `Work` amount. These charges share one thread-local request scope across every decider, including deciders created by recipe and comparison evaluators. Failed reservations latch exhaustion. The top-level publisher checks that latch even if a lower layer catches a refusal. Existing per-decision precision and work limits still apply.

Each request evaluates synchronously on a fresh thread with an 8 MiB stack inside its persistent worker process. The entire numeric scope stays on that thread; nested numeric parallelism is forbidden. Thread teardown discards thread-local intersection caches, so a previous request cannot change expression sharing or charged work. Parallelism is between requests. Wall time includes thread startup and teardown. Only complete external artifacts are cached between requests. This policy bounds accounted work, not CPU instructions or maximum elapsed time. Uninstrumented arithmetic between checkpoints may still be slow. Those are performance bugs to measure, not grounds for a time-dependent geometric answer.

## Driver and evidence

The executable `.claude/skills/verify-catallax/scripts/control-catallax.py` documents commands through `--help` and subcommand `--help`. It supports `doctor`, `run`, `compare`, `report`, `pair`, `legacy` and `cleanup`, each with JSON output. Inputs and evidence are runtime paths. The repository contains no private corpora.

Cache keys include executable SHA-256, source, canonical problem, allowance and request version. Entries contain input, response, manifest, optional STEP and file hashes. Hash validation precedes reuse. Completed entries are checkpoints; an interrupted request is recomputed. Runs preserve input order and duplicate occurrences. `--fresh` bypasses reuse for independence and load tests. Cached wall times are labeled and must not be treated as fresh measurements.

`compare` checks identical ordered inputs, geometric manifests excluding only `stats` and build identity, and exact STEP bytes. Refusals remain visible by stable code. `report` separates outcome counts from slow cases grouped by operation, operand kinds and outcome cause. These groups identify follow-up profiling work; they are not measured call-stack attribution.

`pair` keeps both binaries alive, pins them to the same core, alternates AB and BA, records one warmup then repeated uncached pairs, and retains every response. `legacy` measures the original process-per-question path. Record affinity, build hashes, machine contention and full elapsed run time alongside results. Compare correctness before interpreting timing deltas.
