## 2026-11-07: explicit fillet overflow syntax

Pieces: recipe, features, field, app.

Change: `Operation::Fillet` gains optional `overflow: Option<FilletOverflow>`. The exported `FilletOverflow` accepts `"roll_over"` or `{"keep_neighbours": [FilletSetback, ...]}`. Each setback contains an `EdgeRef`, an endpoint `VertexName` under `at`, and a positive exact `Number` under `setback`.

Why: RECIPE section 4.5.1 already defines station planes, but the parser rejected the policy before either evaluator could report its construction coverage. Version 3 rejects the new field. Versions 4 and 5 validate the station's origin names and positive distance. Endpoint incidence, edge selection and compatible preservation still require geometry.

Migration: Rust constructors supply `overflow: None` to retain their existing version-specific default. Existing serialized recipes are unchanged. Explicit endpoints use the canonical vertex object `{"faces": [...]}`, not the corpus experiment's bare array. Both evaluators refuse `keep_neighbours` with `fillet.unsupported-overflow`; neither substitutes roll-over. This schema change adds no station-plane construction or new checked-export claim.
