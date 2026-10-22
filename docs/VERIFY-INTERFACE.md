# Deterministic verification exchange

`catallax-verify` is an opt-in command. `exact-brep` and `catallax` retain every behavior in [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md), including their wall-clock deadlines. This contract does not replace that exchange.

## Commands and framing

Build with `cargo build --release -p catallax-cli --bins`. `catallax-verify about` prints protocol versions, sources, kernel identity and the default work allowance. `catallax-verify stream` reads one UTF-8 JSON object per line on stdin and flushes one response per line on stdout. EOF shuts down the process. Each line is a complete independent request. Diagnostics go to stderr. A worker processes requests in order; independent workers may process different requests concurrently.

```json
{"format":"catallax-request/1","source":"question","problem":{"format":"exact-brep-question/2","id":"box","operands":{"a":{"shape":"block","size":["1","2","3"]}},"steps":[],"result":"a"},"out":"/tmp/box-proof","max_work":100000000}
```

All fields except `max_work` are required. Unknown envelope fields are invalid. `out` names an owned artifact directory; the command replaces `manifest.json`, `operations.json`, `question.json` and any `answer.step` there. A refusal removes stale STEP. Never assign the same output directory to concurrent requests. `max_work` is an unsigned 64-bit integer; zero is valid and refuses on the first work checkpoint. The default is 100,000,000.

Sources:

| Source | Problem format | Result |
| --- | --- | --- |
| `question` | `exact-brep-question/2` | Full geometric manifest and optional STEP. Parse all legacy fields, but ignore `limits.deadline_ms` during evaluation. |
| `recipe` | `catallax-recipe/3` | Build through `catallax-features`, summarize topology and export STEP with default export options. Recipe properties are not measured by this adapter. |
| `compare` | `catallax-recipe/3` | Run `catallax-compare` with its default options. Report sampled `agree`, `disagree` or `unresolved`, findings and unsupported operations. No STEP output. |
| `field` | `catallax-recipe/3` | Independently compile the field plan. Return candidate source-face count and whether expression identities prove emptiness. No B-rep builder, STEP, sampling or volume measurement. Unsupported operations retain their `field.*` refusal codes. |

A new source implements one evaluation adapter returning a manifest. External corpus adapters export Python `load(Path)` and yield `{"source":...,"problem":...}` records. These adapters are trusted executable code selected explicitly with `--adapter`. Data files cannot select executable adapters.

### Building without measurement

Library callers use `catallax_cli::build_question(&question, &mut decider)` for
question inputs, or `catallax_features::evaluate(&recipe, &mut decider)` for
application recipes. Both return the body and its owning model without computing
proven volume, area or centroid. Neither writes output files.

Request properties later through `catallax_measure::volume`, `area`, `centroid`,
`inertia` or `mass_properties`, passing that model and body with the requested `Accuracy`.
Measurement borrows the model, so a measurement refusal does not discard the
built body or change its construction marks. `build_question` uses the caller's
decider policy, not the question's legacy deadline or output settings.

The question adapter retains its existing measurement defaults. To build through
the verification CLI without measurement, set `"measure":{"properties":[]}` in
the question. Requested measurements use the same enclosure algorithms as the
legacy command. The legacy command and `answer`, `answer_since` and
`answer_deterministic` retain their existing requested-property behavior.

Question sources can explicitly request `"measure":{"properties":["inertia"]}`.
The returned `properties.inertia` is a symmetric 3-by-3 matrix of `{lo,hi}`
enclosures, about the centroid at unit density, in model axes and length units
to the fifth power. Empty bodies return null. Unrequested inertia is omitted.
The default remains volume, area and centroid; existing answers do not compute
second moments or consume an inertia budget. Recipe and comparison outputs are unchanged.

## Response and geometric manifest

A `catallax-response/1` object has `manifest`, `work`, `wall_ms` and the additive `operation_log` field. Consumers must tolerate additional response fields. `work` contains `used`, `limit` and `exhausted`. `wall_ms` is measured elapsed time and never changes an outcome. Work counts are deterministic for the same kernel build, input and allowance. They may change when the implementation changes.

The geometric manifest is `catallax-manifest/3`. For question sources it retains the fields and meanings of manifest/2. Recipe sources return `status`, `exact`, `result`, `marks`, `step`, `refusal` and `kernel`; omitted properties mean unmeasured. A comparison returns `comparison` instead of geometry. `comparison` contains `verdict`, sample counts, ordered findings and `not_comparable`. Findings contain `code`, `detail`, `pass` and diagnostic exact-coordinate `witness`, `brep_faces`, `field_faces` and `gap` strings. These are sampled display checks, never global geometric certificates; `exact` is false for this source. Free-text details and witness formatting are diagnostic.

The surface vocabulary is `plane`, `cylinder`, `elliptic-cylinder`, `cone`, `sphere`, `torus`. The curve vocabulary is `line`, `circle`, `ellipse`, `parabola`, `hyperbola`, `intersection`. Elliptic cylinders use their own tag when their builders land. Parabolas and hyperbolas retain their own tags and finite trimming topology. No old tag is remapped. New analytic kinds require another version change; unknown library kinds must not be published under a convenient existing tag.

`numbers.operation-budget` is a whole-request refusal. It clears partial geometry, marks and STEP output. No deadline refusal originates in this command. Structural input errors have status `invalid-question`; a panic has status `internal-error`. A worker I/O error exits nonzero and is an infrastructure failure. External termination or machine resource exhaustion leaves an incomplete request, never a geometric refusal. A verifier must fail correctness validation for an internal error, missing response or nonzero exit, even if another build fails the same way.

## Operation log

`operation_log` and `operations.json` contain the same `catallax-operation-log/1` object. Its `status`, `refusal` and `marks` describe the final published outcome. `operations` is an execution-ordered array, not a list of planned steps. Failed operations appear once; later operations that never ran do not appear.

| Field | Meaning |
| --- | --- |
| `sequence` | Zero-based position in execution order. |
| `kernel` | `brep` or `field`, the implementation that executed the step. |
| `id` | Original operand or recipe step ID. Publication stages reuse the result ID. |
| `operation` | Recipe operation name, question shape or boolean kind, or `measure`, `summarize`, `step`, `face_measure`. |
| `status` | `answered`, `refused`, or `interrupted` when unwinding prevented completion. |
| `decisions` | Counts from that operation's executing decider, or null for interruption. `total`, `filter`, `refined`, `exact`, `coincident_exact`, `coincident_separation`, `coincident_fact`, `not_exact`, `refused` retain the numeric library's meanings. `max_bits_high_water` is that decider's cumulative precision high-water mark, not a local maximum. |
| `work` | Complete-request counters `before`, `after` and their difference `used`. Observing them charges no work. |
| `marks` | Result marks, including inherited marks, using the manifest mark schema. An empty array is not a local-tolerance certificate. Publication marks also appear at the log root. |
| `local_tolerances` | Null unless this is an opt-in `face_measure` operation. That operation carries an array containing its face's definition-deviation certificate, as described below. Marks and STEP deviation are not relabelled as local tolerances. |
| `refusal` | Null or an object with stable `code` and diagnostic `detail`. |
| `wall_ms` | Observed operation time. Exclude it when comparing deterministic logs. |

Instrumentation covers each executed B-rep recipe step, each field compilation step, question operands and booleans, and B-rep publication stages. Comparison sources include the build steps they execute, but comparison sampling remains accounted in the response's whole-request `work` and `wall_ms`, not attributed to a recipe step. Validation and internal deciders can also charge work outside the listed spans. Do not sum span counters as a substitute for whole-request accounting. A comparison's existing unsupported-operation preflight may return before executing any step.

The log records decision settlement counts, not expressions or proof certificates for individual numeric calls. It does not recompute decisions to obtain diagnostics. Recording is disabled outside the verification exchange and adds no work charges. No numeric precision, geometry, outcome or STEP export option changes.

A budget refusal discards published geometry and marks as before. The log keeps its executed prefix for diagnosis, with the final `numbers.operation-budget` refusal at its root. Prefix marks do not authorize publishing a partial body. A panic remains `internal-error` even when interrupted operations were recorded.

## Whole-face local tolerances

The `recipe` adapter also accepts a `catallax-face-check/1` problem. It wraps an
unchanged recipe with typed per-face definitions and an optional scalar
parameter edit. Requests without this wrapper keep their existing fields,
work, geometry and STEP bytes. Neither the legacy exchange nor recipe/3 changes.

This portable question exercises a cylinder side as the rolling surface between
two perpendicular planes whose outward unit offsets meet on the z axis:

```json
{
  "format": "catallax-question/1",
  "source": "recipe",
  "max_work": 100000000,
  "problem": {
    "format": "catallax-face-check/1",
    "recipe": {
      "format": "catallax-recipe/3", "name": "rolling-cylinder",
      "steps": [{
        "id": "tube", "op": "cylinder", "radius": "1", "height": "3",
        "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}
      }],
      "result": "tube"
    },
    "definitions": [{
      "face": {"step": "tube", "role": {"side": {"segment": "wall"}}},
      "definition": {
        "kind": "fillet", "radius": "1",
        "path": {"kind": "line", "origin": ["0", "0", "0"], "direction": ["0", "0", "1"]}
      }
    }],
    "options": {"width": "0.0000001", "max_cells": 4096, "path_segments": 64},
    "edit": {"step": "tube", "parameter": "radius", "value": "1.1"}
  }
}
```

Save it outside the repository and use `catallax-verify replay FILE --out DIR`
or the driver's `run` command. The initial face has zero deviation, not a null
tolerance. The edit rebuilds the recipe and resolves the same origin references
on the edited result. It reports names after the edit, not a second measurement
with an implicitly changed definition.

`definitions` must be nonempty. Each `face` is a recipe `OriginPattern`; every
matching result face receives that definition. Empty selections and duplicate
definitions of one face refuse with `cli.face-definition`. Exact faces need no
record. Definitions are:

| `kind` | Fields and checked definition |
| --- | --- |
| `fillet` | `path`, positive `radius`. Check `abs(distance(point, path) - radius)` over the whole face. The path is the ball-centre path, not either original rolling support. |
| `offset` | `original`, non-negative `distance`. Check `abs(distance(point, original) - distance)` against the original analytic support, not the whole solid's boundary or Boolean residual. |

Path records use these exact fields:

| `path.kind` | Fields |
| --- | --- |
| `line` | `origin`, nonzero `direction`. The whole line. |
| `circle` | `center`, `normal`, positive `radius`, perpendicular `reference`. The whole circle. |
| `ellipse` | `center`, `normal`, perpendicular `major_direction`, `major`, `minor`, with `0 < minor < major`. The whole ellipse. |
| `edge` | `between`, two origin patterns resolving exactly one edge in the result. Uses its exact carrier, selected branch and arc endpoints. |
| `intersection` | `surfaces`, two exact supporting-surface records already at the outward offset; `component`, the zero-based curve component in deterministic intersection order; optional exact `from` and `to` points. Absent endpoints select the whole closed curve, or the complete line for a line component. Non-closed arcs otherwise require their endpoints. |

Supporting-surface records use `kind` with these fields:

| `original.kind` or `surfaces[].kind` | Fields |
| --- | --- |
| `plane` | `origin`, `normal`, perpendicular `reference` |
| `cylinder` | `origin`, `axis`, positive `radius`, perpendicular `reference` |
| `cone` | `apex`, `axis`, positive `slope`, perpendicular `reference` |
| `sphere` | `center`, positive `radius`, `axis`, perpendicular `reference` |
| `torus` | `center`, `axis`, `major`, `minor`, perpendicular `reference`, with `0 < minor < major` |
| `elliptic_cylinder` | `origin`, `axis`, perpendicular `major_direction`, `major`, `minor`, with `0 < minor < major` |

Points and directions are arrays of three recipe-number strings. All dimensions
and controls use exact recipe numbers; limits are positive integer counts.
Unknown fields refuse. The wrapper is available only through `source=recipe`;
the B-rep builder still refuses unimplemented fillet and offset operations.
Declaring a definition measures a constructed face and does not construct a
new feature.

The successful manifest adds `face_measurements`, format
`catallax-face-measurements/1`:

- `faces` lists each measured topology-local `face`, its complete resolved
  `origins`, definition `kind`, `deviation_mm.lower` and `deviation_mm.upper`,
  `local_tolerance_mm`, `converged`, covering `cells`, exact `witnesses` and
  diagnostic `coverage`. Bound endpoints and the local tolerance are exact
  rational strings in millimetres. `local_tolerance_mm` is the upper endpoint.
- `display` reports observed `first_mesh_ms` from request recipe processing
  through the first complete B-rep display mesh, plus `vertices`, `triangles`
  and `chord_tolerance_mm`. It excludes definition checking, the edit rebuild,
  summary, STEP export, transport and UI presentation. Display chord tolerance
  is not the definition tolerance.
- `edit` is null or contains `step`, `parameter`, canonical `value` and a
  `resolved` entry per declared reference. Each entry has `reference`, `status`
  of `resolved` or `missing`, and all matching edited `faces` with their
  complete origins. The edit must name an existing scalar step parameter.
  Geometry is rebuilt through the normal recipe evaluator.

`options` is optional and defaults to the values in the example. `width` is the
requested width of the worst-deviation enclosure, not a geometric acceptance
threshold. A cell limit returns a valid bound with `converged=false`. Upper
certificates cover complete supporting chart rectangles, which may include
trimmed-away regions. Lower bounds use exact witnesses on the actual trimmed
face. See [FIELD.md](FIELD.md#whole-face-definition-measurements) for the
derivative-cover and closest-path proofs and propagated refusal limits.

Each measured face adds an executing field `face_measure` operation whose
`local_tolerances` contains that face's certificate. Unrelated operations retain
null. The measurement is not a tolerance-core cap or merge decision, and it does
not turn sampled cross-kernel agreement into a whole-part distance certificate.

`first_mesh_ms` is an observation. Exclude it, operation `wall_ms`, manifest
`stats` and build identity when comparing deterministic runs. Retain every
geometric certificate, origin, work count and strict STEP byte.

## Application-exported questions

An application writes a UTF-8 JSON object with exactly these four fields:

```json
{"format":"catallax-question/1","source":"recipe","problem":{"format":"catallax-recipe/3","name":"exported-block","steps":[{"id":"base","op":"block","size":["4","5","6"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"base"},"max_work":12345678}
```

`source` is one of the exchange adapters above. `problem` is the complete original input, with its original IDs and selected result. `max_work` is the original unsigned 64-bit allowance and is required, including when zero. The format contains no machine-specific output path, expected answer, timestamp or executable selection.

`catallax-verify replay QUESTION.json --out DIR` reads this format, binds the artifact directory and executes the same adapter as `stream`. Every completed exchange also writes its portable question to `question.json`. The batch driver's `run` and `pair` read portable questions and preserve their recorded allowance instead of substituting the command's default. A failing feature can be exported as a prefix only if the application originally ran that prefix; changing the evaluated recipe can change work and outcomes.

Replay equality means the same published geometric manifest, strict STEP bytes and complete-request work count on the same kernel revision and allowance. Wall times are observations, not part of equality. Application command history, screenshots and UI state remain the application's responsibility.


## Work policy

Each `Decider` and `Work` cancellation checkpoint costs one unit. Accepted arithmetic reservations also cost their existing deterministic `Work` amount. These charges share one thread-local request scope across every decider, including deciders created by recipe and comparison evaluators. Failed reservations latch exhaustion. The top-level publisher checks that latch even if a lower layer catches a refusal. Existing per-decision precision and work limits still apply.

Each request evaluates synchronously on a fresh thread with an 8 MiB stack inside its persistent worker process. The entire numeric scope stays on that thread; nested numeric parallelism is forbidden. Thread teardown discards thread-local intersection caches, so a previous request cannot change expression sharing or charged work. Parallelism is between requests. Wall time includes thread startup and teardown. Only complete external artifacts are cached between requests. This policy bounds accounted work, not CPU instructions or maximum elapsed time. Uninstrumented arithmetic between checkpoints may still be slow. Those are performance bugs to measure, not grounds for a time-dependent geometric answer.

## Driver and evidence

The executable `.claude/skills/verify-catallax/scripts/control-catallax.py` documents commands through `--help` and subcommand `--help`. It supports `doctor`, `run`, `compare`, `report`, `pair`, `legacy`, `parts` and `cleanup`, each with JSON output. Inputs and evidence are runtime paths. The repository contains no private corpora.

Cache keys include executable SHA-256, source, canonical problem, allowance and request version. Entries contain input, response, manifest, optional STEP and file hashes. Hash validation precedes reuse. Completed entries are checkpoints; an interrupted request is recomputed. Runs preserve input order and duplicate occurrences. `--fresh` bypasses reuse for independence and load tests. Cached wall times are labeled and must not be treated as fresh measurements.

`compare` checks identical ordered inputs, geometric manifests excluding `stats`, build identity and opt-in `face_measurements.display.first_mesh_ms`, and exact STEP bytes. Refusals remain visible by stable code. `report` separates outcome counts from slow cases grouped by operation, operand kinds and outcome cause. These groups identify follow-up profiling work; they are not measured call-stack attribution.

`pair` keeps both binaries alive, pins them to the same core, alternates AB and BA, records one warmup then repeated uncached pairs, and retains every response. `legacy` measures the original process-per-question path. Record affinity, build hashes, machine contention and full elapsed run time alongside results. Compare correctness before interpreting timing deltas.

## Target-part scoreboard

```sh
python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py parts parts \
  --binary target/release/catallax-verify --out /tmp/catallax-parts \
  --max-work 1000000 --core 2 --occt-python /tmp/occt/bin/python --json
```

Use an allowed core from `doctor`, or omit `--core`. Install `cadquery-ocp` in a scratch virtual environment outside the repository for `--occt-python`. A missing binding is `step.binding-unavailable`, not a passing check. OpenCASCADE never runs inside the kernel or becomes a runtime dependency of Catallax.

`parts` reads every recipe recursively in filename order and writes `parts.json` with format `catallax-parts/1`. Each row retains the relative `recipe_file`, SHA-256, name, and complete independent `recipe`, `field` and `compare` responses, including work, wall time, refusal codes and operation logs. All three adapters run even if one refuses. The report is checkpointed after each completed row; its request, response, manifest, log, portable question and any STEP remain under `cases/`. Runs are fresh, not cached.

For each exported STEP, a separate process opens the file in OpenCASCADE, transfers its roots, checks validity without repair, compares unique face count and records imported volume. The STEP hash and binding version accompany the result. `step_check.status=passed` covers import, topology and face count only. The recipe adapter does not measure kernel volume; `volume_comparison` explicitly reports `step.kernel-volume-unmeasured`. No volume-deviation or full target-part acceptance claim follows from a passing import. No STEP is `step.no-export` with its build refusal, not a silently skipped part.

The report labels machine load and affinity. Shared-machine observations are not quiet-core speed evidence. `application_check` remains `parts.application-loop-required`; the scoreboard does not pretend to have edited a part in SolidDraft. Sampled comparison still is not a global geometric certificate. Nonzero exit reports invalid input, internal error, proven comparison disagreement or a failed STEP import; ordinary coverage refusals and unresolved comparisons stay visible without making report generation fail.

`.github/workflows/parts.yml` runs this report on `main` and relevant pull requests at a fixed 1,000,000-unit allowance and uploads the full bundle, even on failure. It is a diagnostic workflow separate from the full-suite `check` gate. Checkpoint reports use the same explicit allowance so a coverage refusal is never mistaken for impossible geometry.
