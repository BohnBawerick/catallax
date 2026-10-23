# Deterministic answers

Build exact geometry or compare the two recipe evaluators without a wall-clock cutoff choosing the result.

## Sub-features

- `answer.question`: question/2 produces manifest/3 and requested STEP.
- `answer.recipe`: recipe/3 produces topology and STEP.
- `answer.compare`: recipe/3 produces a sampled comparison report.
- `answer.field-display`: recipe/3 produces the field's own `display.obj`, a proven two-sided boundary distance bound and `face_measurements.display.first_mesh_ms`. Unsupported features and charts refuse without B-rep fallback.
- `answer.budget`: zero work refuses and leaves no STEP.
- `answer.log`: each completed response records executed operations with decision settlement counts, marks, work, refusal codes and diagnostic timings.
- `answer.export`: portable question/1 replays the original source and work allowance.
- `answer.parts`: P0–P5 builds, independent field compilation, sampled comparison and STEP import report.
- `answer.face-definitions`: opt-in recipe face checks report whole-face two-sided deviation, local tolerances, first-mesh timing and origins after a scalar edit.

## How to get to it (user POV)

Run the driver's `run` command with a question, recipe or directory. Select `--source compare` for the two-kernel path. Direct streaming uses `catallax-verify stream` with the documented envelope.

## Driving it with control-catallax

Preconditions: release build, doctor passed, `PROOF` names an external directory.

- Question: run `control-catallax.py run .claude/skills/verify-catallax/references/examples/block.json --binary target/release/catallax-verify --out "$PROOF/question" --json`. Require one answered result and an `answer.step` artifact despite the example's zero legacy deadline.
- Recipe: run `control-catallax.py run .claude/skills/verify-catallax/references/examples/recipe.json --source recipe --binary target/release/catallax-verify --out "$PROOF/recipe" --json`. Require answered topology and STEP.
- Comparison: repeat the recipe command with `--source compare` and a separate output directory. Inspect `manifest.comparison.verdict` and findings. Require no STEP.
- Field display: repeat the recipe command with `--source field` and a separate output directory. Require `display.obj`, actual vertex and triangle counts, `face_measurements.display.first_mesh_ms` and exact rational `deviation_bound_mm`. The empty `face_measurements.faces` list is not a fillet-definition measurement. Replay `parts/P3-filleted-boss.json` with this source and require `field.unsupported-operation` until the field evaluator implements fillets, with neither stale OBJ nor STEP. Compare independent one-, two- and four-worker runs, retaining strict OBJ bytes as well as certificates.
- Budget: repeat the question command with `--max-work 0` and a separate directory. Require `numbers.operation-budget` and no STEP.
- Export: locate the recipe run's `question.json`, then run `target/release/catallax-verify replay QUESTION --out "$PROOF/replayed"`. Require the same manifest, work and STEP bytes. Compare executed operation entries after removing only `wall_ms`.
- Scoreboard: run `control-catallax.py parts parts --binary target/release/catallax-verify --out "$PROOF/parts" --max-work 1000000 --occt-python /tmp/occt/bin/python --json`. Require every recipe in `parts.json`, including typed refusals from both implementations. Check `step_check` for actual OpenCASCADE import status and coverage limits.
- Face definitions: save the portable cylinder example from `docs/VERIFY-INTERFACE.md` under an external proof directory and run it through `control-catallax.py run FILE --binary target/release/catallax-verify --out "$PROOF/faces" --json`. Require zero lower and upper deviation for the declared rolling cylinder, a non-null `face_measure.local_tolerances` entry, actual first-mesh timing, resolved wall origins after the radius edit and an `answer.step` artifact. For refinement, translate a torus relative to its declared circle path and increase `options.max_cells`; require an enclosing interval that tightens, not a sample maximum. Deterministic comparisons exclude only the diagnostic first-mesh clock in addition to existing stats and build identity.

## Gotchas

Use the full driver path from SKILL.md if it is not on PATH. Comparison agreement is sampled. A valid unsupported recipe can honestly refuse. Recipe export uses default export options and does not measure mass properties.
