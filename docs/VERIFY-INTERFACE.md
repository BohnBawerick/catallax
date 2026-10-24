# Deterministic verification exchange

`catallax-verify` is an opt-in command. `exact-brep` and `catallax` retain every behavior in [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md), including their wall-clock deadlines. This contract does not replace that exchange.

## Commands and framing

Build with `cargo build --release -p catallax-cli --bins`. `catallax-verify about` prints protocol versions, sources, kernel identity and the default work allowance. `catallax-verify stream` reads one UTF-8 JSON object per line on stdin and flushes one response per line on stdout. EOF shuts down the process. Each line is a complete independent request. Diagnostics go to stderr. A worker processes requests in order; independent workers may process different requests concurrently.

```json
{"format":"catallax-request/1","source":"question","problem":{"format":"exact-brep-question/2","id":"box","operands":{"a":{"shape":"block","size":["1","2","3"]}},"steps":[],"result":"a"},"out":"/tmp/box-proof","max_work":100000000}
```

All fields except `max_work` are required. Unknown envelope fields are invalid. `out` names an owned artifact directory; the command replaces `manifest.json`, `operations.json`, `question.json`, `answer.step` and `display.obj` there. A refusal removes stale geometry artifacts. Never assign the same output directory to concurrent requests. `max_work` is an unsigned 64-bit integer; zero is valid and refuses on the first work checkpoint. The default is 100,000,000.

Sources:

| Source | Problem format | Result |
| --- | --- | --- |
| `question` | `exact-brep-question/2` | Full geometric manifest and optional STEP. Parse all legacy fields, but ignore `limits.deadline_ms` during evaluation. |
| `recipe` | `catallax-recipe/3`, `catallax-recipe/4` | Build through `catallax-features`, summarize topology and export STEP with default export options. Recipe properties are not measured by this adapter. |
| `compare` | `catallax-recipe/3`, `catallax-recipe/4` | Run `catallax-compare` with its default options. Report sampled `agree`, `disagree` or `unresolved`, findings and unsupported operations. No STEP output. |
| `field` | `catallax-recipe/3`, `catallax-recipe/4` | Independently compile the field plan and build its native bounded display mesh in `display.obj`. Report candidate source-face count, proven emptiness and `face_measurements.display`. No B-rep builder, STEP or volume measurement. Unsupported operations and geometry retain their specific refusal codes, with no display fallback. |
| `name-edit` | `catallax-name-edit/1` | Rebuild `recipe` and `edited_recipe`, resolve every original result-face origin and physical-edge name, and check the independently compiled field registries. Report `name_edit`; no STEP or display. |

A new source implements one evaluation adapter returning a manifest. External corpus adapters export Python `load(Path)` and yield `{"source":...,"problem":...}` records. These adapters are trusted executable code selected explicitly with `--adapter`. Data files cannot select executable adapters.

### Parameter-edit names

`catallax-name-edit/1` contains `format`, `recipe` and `edited_recipe`, both complete recipes. Result and ordered step IDs and operation kinds must stay fixed. An unchanged input is `names.no-edit`. All four builds share the request's unchanged deterministic allowance.

`name_edit.status` is `passed`, `failed`, `waiting` or `invalid`. Face references resolve by origin containment, including split pieces and merged origin sets. Edge references resolve against both incident faces' complete origin sets, not arena handles or only their current first origins. Missing references and unnamed physical entities fail the check. The report counts checked faces and edges and lists missing names. Field coverage is explicit: result-face origins resolve in each independent definition registry, not a proof of surviving boundary domains. Sampled geometry comparison remains a separate check.

A build refusal reports `waiting` with kernel, before/after phase, refusal code and detail. B-rep refusals include the exact step ID and operation, so a shell or chamfer still being built cannot become a passing edit check. Invalid envelopes and structural recipes report `invalid`.

### Field display

`field` runs the field's own `catallax_field::display::mesh_plan` path. The CLI
does not dispatch on feature kinds or obtain a B-rep answer. A future field
feature supplies its definition and certified native display path through the
field library; the verification adapter does not change.

An answered request writes `display.obj`, including an empty mesh for a proven
empty result. Its manifest adds `face_measurements`, format
`catallax-face-measurements/1`, with empty `faces`, null `edit` and `display`:

| Display field | Meaning |
| --- | --- |
| `first_mesh_ms` | Elapsed request recipe processing, field compilation and native meshing, through the first complete display mesh. Same start and completion boundary as the B-rep face-check clock. |
| `file` | `display.obj`, the actual binary64 display geometry. |
| `vertices`, `triangles` | Actual mesh counts. |
| `chord_tolerance_mm` | Requested display tolerance, currently 0.01 mm. |
| `deviation_bound_mm` | Proven two-sided upper distance between the mesh and the field-defined result boundary, including trim coverage and coordinate rounding. Exact rational string in millimetres. |
| `coverage` | Diagnostic description of that certificate. |
| `cells` | Deterministic visited parameter-cell count. |

The clock excludes writing OBJ, response transport and UI presentation, just
as the B-rep clock excludes summary, export and presentation. This is not
time to a rendered application frame. The display bound is not a whole-face
fillet construction measurement; empty `faces` does not certify a fillet.
`exact` describes the field definition, not the triangle approximation.

Meshing adds an executed field `display` operation. Invalid controls refuse
with `field.display-options`; unsupported native charts or exact frames use
`field.display-unsupported`; missing compact bounds use
`field.display-unbounded`; resource/refinement limits use
`field.display-budget`; coordinate range or rounding failures use
`field.display-precision`. Numeric refusals propagate. Every refusal removes
`display.obj` and `answer.step`. Elementary fillets compile independently;
unsupported contact domains retain their specific construction refusals.
No partial mesh, B-rep fallback or sampled deviation estimate is published.

Executed mesh attempts also publish `display_diagnostics`, format
`catallax-display-diagnostics/1`, even on a shared numeric-budget refusal.
This is an observed prefix, never partial geometry or a distance certificate.
`total` and `common` contain `cells`, `max_depth`, `decisions` and `work`.
`work` is accepted complete-request work; reading a counter charges nothing.
`common` covers preparation and shared indexing not assigned to a support.
Decisions count the executing caller's decider, as in the operation log.

`supports` groups coincident source faces. Each entry has registry `faces`,
`origins`, `total`, and included `classification`, `queries` and `emission`
stage costs. Count shared support work once. `faces` records each source
registry handle and origin, shared `candidate_cells` and `max_depth`, and
its executed boundary `queries` cost. A coactive name may have no separate
query when another name already proved the same patch. Face query costs are
included in support query costs, not added to them. All these observations
are deterministic and remain part of independent-worker comparisons.


The driver accepts `--source field` for recipe inputs. `compare` excludes
`first_mesh_ms` but compares display certificates, work for identical builds,
strict `display.obj` bytes and strict STEP bytes.


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

`comparison.sampled_distance_upper_mm` is the maximum certified distance upper
bound over checked samples, including forward mesh error. It is an exact rational
string in millimetres, or `null` if any required sample distance is unbounded.
It is not a whole-part Hausdorff bound. Trimmed-face queries can establish finite
upper bounds on curved Boolean and swept faces without changing the comparison's
accuracy or promoting residual values to true distance.

The surface vocabulary is `plane`, `cylinder`, `elliptic-cylinder`, `cone`, `sphere`, `torus`. The curve vocabulary is `line`, `circle`, `ellipse`, `parabola`, `hyperbola`, `intersection`. Elliptic cylinders use their own tag when their builders land. Parabolas and hyperbolas retain their own tags and finite trimming topology. No old tag is remapped. New analytic kinds require another version change; unknown library kinds must not be published under a convenient existing tag.

`numbers.operation-budget` is a whole-request refusal. It clears partial geometry, marks and STEP output. No deadline refusal originates in this command. Structural input errors have status `invalid-question`; a panic has status `internal-error`. A worker I/O error exits nonzero and is an infrastructure failure. External termination or machine resource exhaustion leaves an incomplete request, never a geometric refusal. A verifier must fail correctness validation for an internal error, missing response or nonzero exit, even if another build fails the same way.

## Operation log

`operation_log` and `operations.json` contain the same `catallax-operation-log/1` object. Its `status`, `refusal` and `marks` describe the final published outcome. `operations` is an execution-ordered array, not a list of planned steps. Failed operations appear once; later operations that never ran do not appear.

| Field | Meaning |
| --- | --- |
| `sequence` | Zero-based position in execution order. |
| `kernel` | `brep` or `field`, the implementation that executed the step. |
| `id` | Original operand or recipe step ID. Publication stages reuse the result ID. |
| `operation` | Recipe operation name, question shape or boolean kind, or `measure`, `summarize`, `step`, `face_measure`, `display`. |
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

## Certified recipe bodies (opt-in)

An existing `catallax-request/1` with `source:"recipe"` can carry a
`catallax-certified-recipe/1` problem. The nested `recipe` remains unchanged
`catallax-recipe/3`; ordinary recipe, question, comparison and face-check inputs
retain manifest/3 and their existing semantics. `about` advertises the additive
`certified_recipe` and `certified_manifest` versions without changing its default.

The public [certified recipe example](../.claude/skills/verify-catallax/references/examples/certified-recipe.json)
constructs the approved oblique-circle/clip source and its inward elliptic
normal-offset prism. Run the raw example through `control-catallax.py run`;
the adapter selects `recipe`. The wrapper contains:

| Field | Contract |
| --- | --- |
| `recipe` | Required recipe/3 object. |
| `outputs.mesh` | Optional `{chord_tolerance,max_vertices}`; writes `display.obj`. Chord tolerance is in mm and concerns the represented model, not its exact target. |
| `outputs.step` | Optional `{file_tolerance,decimal_places}`; writes `answer.step` only after certification. The nonnegative tolerance bounds file-to-model distance in mm. |
| `properties` | Optional `{relative,absolute,inertia}`; nonnegative requested enclosure widths. `inertia` defaults false. |
| `field_points` | Optional array of exact three-component points in mm. An empty array requests no field output. |
| `edit` | Optional `{step,parameter,value}` scalar edit, using the existing face-check edit contract. Rebuilds a separate snapshot; does not mutate the initial body or its artifacts. |

`outputs` itself is required and may be `{}`. Numeric inputs use recipe `Number`
syntax, including exact rational strings. Unknown wrapper/option fields are
invalid input.

Only this opt-in problem produces `catallax-manifest/4`. Its `result` uses the
existing analytic kind tags plus genuine `bspline` surface and curve tags.
The first supported spline body has one side, two complete planar cap disks,
two shared physical ring edges and no model vertices. Its finite spline domain
and knots are reported in `certificate.spline`; `u=0` and `u=4` identify a C1
closed boundary, not an extra physical model edge. Each cap and side region
records its actual face, origins, complete trim domain and evidence binding.

The complete `certificate` reports revision, representation, target and scope
identities; exact regions and local exceptions; the evidence dependency graph;
proof kinds; the decide-once branch/evidence records; and any proved uniform
signed-distance bound. Exceptions retain their measured enclosure, source,
hard limit and proved neighbouring edge/face size lower bounds. Settings record
the separate 1 nm modelling tolerance, 100 nm engineering target, 1 micrometre
kernel cap and 10 micrometre import cap. A missing uniform metric is null, not
zero. `exact` requires an entirely exact region cover and no bounded decisions;
an empty marks array does not establish it. These wire records report checked
claims: hashes and parsed JSON alone never confer geometric authority. Replay
the retained recipe and recheck geometry to reproduce that authority.

Requested consumers return separate `display`, `properties`, `field`, `step`
and `edit` objects with `status:"answered"` or `"refused"`; unrequested consumers
are null. A consumer refusal preserves the complete valid body and names. A
whole-request `numbers.operation-budget` refusal instead retains manifest/4
but clears geometry, certificate, consumer results and output artifacts. Panics
and I/O failures retain the existing infrastructure/error rules.

Display face groups and STEP `file_region` names are `face_<snapshot handle>`.
Their region records distinguish represented-model chord/file errors from
composed output-to-original-definition bounds. File accuracy never enlarges a
kernel tolerance. STEP uses an actual degree-(3,1) B-spline surface and its
shared rounded boundary curves; export-only vertices/edges are reported
separately from model topology. Unsupported certified export coverage refuses.

Properties and field queries concern the retained **original definition**, not
the fitted carrier. Properties retain outward decimal enclosures and must meet
their requested width without adding the carrier error. Field results report
true signed distance, location, unsigned distance and normal/ambiguity.
Certificate, field and output-error endpoints are exact rational strings;
unbounded field endpoints are `"infinity"` or `"-infinity"`, never JSON null.
The property DTO retains its existing outward-decimal representation and units.

Edited results include their new certificate, names and origin/edge-name matches
for old references. Handles are snapshot-local; inspect all matches rather than
assuming a surviving reference means an unchanged handle or a unique match.
`replay:"question.json"` names the portable request including its work allowance.
For deterministic comparison, exclude diagnostic clocks and build identity as
before, but compare certificates, knots, domains, consumer outcomes, names and
artifact bytes.

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
the B-rep builder supports elementary cylinder and torus fillet strips but still
refuses unimplemented general contacts and offsets. Declaring a measurement
definition does not construct a new feature.

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

`pair` keeps both binaries alive, pins them to the same core, alternates AB and BA, records one warmup then repeated uncached pairs, and retains every response. Each pair records `load_average_before` and `load_average_after`; these are the operating system's one-, five- and fifteen-minute averages, not per-process CPU measurements. `legacy` measures the original process-per-question path. Record affinity, build hashes, machine contention and full elapsed run time alongside results. Compare correctness before interpreting timing deltas.

## Target-part scoreboard

```sh
python3 -B .claude/skills/verify-catallax/scripts/control-catallax.py parts parts \
  --binary target/release/catallax-verify --out /tmp/catallax-parts \
  --max-work 1000000 --core 2 --occt-python /tmp/occt/bin/python --json
```

Use an allowed core from `doctor`, or omit `--core`. Install `cadquery-ocp` in a scratch virtual environment outside the repository for `--occt-python`. A missing binding is `step.binding-unavailable`, not a passing check. OpenCASCADE never runs inside the kernel or becomes a runtime dependency of Catallax.

`parts` reads every recipe recursively in filename order and writes `parts.json` with format `catallax-parts/1`. Each row retains the relative `recipe_file`, SHA-256, name, and complete independent `recipe`, `field` and `compare` responses, including work, wall time, refusal codes and operation logs. All three adapters run even if one refuses. The report is checkpointed after each completed row; its request, response, manifest, log, portable question and any STEP remain under `cases/`. Runs are fresh, not cached.

`parts/edit-checks.json`, format `catallax-part-edits/1`, is a specification rather than a recipe. Its `cases` map relative recipe filenames to edits with `step`, a nonempty `path` of object keys or array indices inside that step, and an exact numeric string `value`. The report rebuilds the complete edited recipe through `name-edit`, saves its full exchange artifacts under `edit/`, and adds `edit_check` to the row. P1 changes boss radius 8 to 9; P2 changes chamfer setback 1 to 0.75; P4 changes shell thickness 2 to 1.5. Other rows report `names.no-edit-specified`. `edit_failures` counts failed or invalid checks and makes the command exit nonzero. Waiting checks stay visible without claiming the roadmap exit.

Configured name-edit checks run and checkpoint before display or sampled comparison. The top-level `edit_checks` map retains each completed check even if an unrelated geometry row is still running. Report `status` is `running`, `completed` or `interrupted`; partial reports do not imply completed geometry coverage. When request exhaustion removes the nested naming report, the checkpoint retains the refusal and first `blocked_operation` from the operation log. The allowance is unchanged.

For each exported STEP, a separate process opens the file in OpenCASCADE, transfers its roots, checks validity without repair, compares unique face count and records imported volume. The STEP hash and binding version accompany the result. `step_check.status=passed` covers import, topology and face count only. The recipe adapter does not measure kernel volume; `volume_comparison` explicitly reports `step.kernel-volume-unmeasured`. No volume-deviation or full target-part acceptance claim follows from a passing import. No STEP is `step.no-export` with its build refusal, not a silently skipped part.

The report labels machine load and affinity. Shared-machine observations are not quiet-core speed evidence. `application_check` remains `parts.application-loop-required`; the scoreboard does not pretend to have edited a part in SolidDraft. Sampled comparison still is not a global geometric certificate. Nonzero exit reports invalid input, internal error, proven comparison disagreement, a failed STEP import or a failed/invalid name edit; ordinary coverage refusals, waiting edits and unresolved comparisons stay visible without making report generation fail.

`.github/workflows/parts.yml` runs this report on `main` and relevant pull requests at a fixed 1,000,000-unit allowance and uploads the full bundle, even on failure. It is a diagnostic workflow separate from the full-suite `check` gate. Checkpoint reports use the same explicit allowance so a coverage refusal is never mistaken for impossible geometry.

### Released-tool head-to-head

`control-catallax.py head-to-head` adds a development-only `catallax-head-to-head/1` report beside the existing scoreboard. It uses the `catallax-cli` example `head_to_head`, not an application or kernel runtime entry with added measurement. The [rerun guide](../.claude/skills/verify-catallax/scripts/head-to-head/README.md) defines installation, pinned reference release, authored P0–P5 scripts, matched pairs, phase scopes and independent property and STEP checks.

The raw report retains construction, complete recipe edit, display, requested measurement and export phases, refusal timing or its explicitly scoped upper bound, work and operations rebuilt. Missing reference work counters remain null. It preserves incorrect answers, failed commands, incomplete requests, missing names and unavailable checks rather than omitting targets. OpenCASCADE independently constructs the reference geometry; its numerical comparisons and sampled membership checks are diagnostics, not kernel certificates. A completed report is not a claim that every target built or met its click budget. `correctness_failures` makes detected property, STEP, membership and name failures return a nonzero command status while retaining the report.

