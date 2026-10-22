# Deterministic answers

Build exact geometry or compare the two recipe evaluators without a wall-clock cutoff choosing the result.

## Sub-features

- `answer.question`: question/2 produces manifest/3 and requested STEP.
- `answer.recipe`: recipe/3 produces topology and STEP.
- `answer.compare`: recipe/3 produces a sampled comparison report.
- `answer.budget`: zero work refuses and leaves no STEP.
- `answer.log`: each completed response records executed operations with decision settlement counts, marks, work, refusal codes and diagnostic timings.
- `answer.export`: portable question/1 replays the original source and work allowance.
- `answer.parts`: P0–P5 builds, independent field compilation, sampled comparison and STEP import report.

## How to get to it (user POV)

Run the driver's `run` command with a question, recipe or directory. Select `--source compare` for the two-kernel path. Direct streaming uses `catallax-verify stream` with the documented envelope.

## Driving it with control-catallax

Preconditions: release build, doctor passed, `PROOF` names an external directory.

- Question: run `control-catallax.py run .claude/skills/verify-catallax/references/examples/block.json --binary target/release/catallax-verify --out "$PROOF/question" --json`. Require one answered result and an `answer.step` artifact despite the example's zero legacy deadline.
- Recipe: run `control-catallax.py run .claude/skills/verify-catallax/references/examples/recipe.json --source recipe --binary target/release/catallax-verify --out "$PROOF/recipe" --json`. Require answered topology and STEP.
- Comparison: repeat the recipe command with `--source compare` and a separate output directory. Inspect `manifest.comparison.verdict` and findings. Require no STEP.
- Budget: repeat the question command with `--max-work 0` and a separate directory. Require `numbers.operation-budget` and no STEP.
- Export: locate the recipe run's `question.json`, then run `target/release/catallax-verify replay QUESTION --out "$PROOF/replayed"`. Require the same manifest, work and STEP bytes. Compare executed operation entries after removing only `wall_ms`.
- Scoreboard: run `control-catallax.py parts parts --binary target/release/catallax-verify --out "$PROOF/parts" --max-work 1000000 --occt-python /tmp/occt/bin/python --json`. Require every recipe in `parts.json`, including typed refusals from both implementations. Check `step_check` for actual OpenCASCADE import status and coverage limits.

## Gotchas

Use the full driver path from SKILL.md if it is not on PATH. Comparison agreement is sampled. A valid unsupported recipe can honestly refuse. Recipe export uses default export options and does not measure mass properties.
