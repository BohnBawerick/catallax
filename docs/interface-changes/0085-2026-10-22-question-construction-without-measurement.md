## 2026-10-22: question construction without measurement

Pieces: CLI library, application callers, measurement
Change: `catallax_cli::build_question(&Question, &mut Decider)` returns a `BuiltQuestion` with its owning `Model`, result `body`, construction `marks` and boolean `per_step` statistics. It performs no property integration, topology summary or export. `BuildError` retains the refusal and optional boolean step and stage.
Why: an application needs a usable body before requesting proven volume, area or centroid. A later measurement refusal must not discard that body.
Migration: library consumers that need a body use `build_question`, then call `catallax_measure::volume`, `area` or `mass_properties` on demand. Recipe applications keep using the already unmeasured `catallax_features::evaluate` path. The legacy commands, answer functions, verification defaults, requested enclosures and STEP contract are unchanged. Verification questions can continue opting out with an empty `measure.properties` list.

