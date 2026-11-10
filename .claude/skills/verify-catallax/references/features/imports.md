# Imported STEP fields

Build the release CLI and pass doctor as in SKILL.md. Use an external proof directory and the original work allowance.

## Input

Use a `catallax-import-recipe/1` wrapper with exactly `format`, `recipe` and `files`. `recipe` is an ordinary recipe with `import_step` operations. `files` maps each saved SHA-256 to complete UTF-8 STEP text. The CLI checks every supplied hash. It opens no source paths from recipe file labels. See [the exchange contract](../../../../../docs/VERIFY-INTERFACE.md).

## Drive

Run the wrapper through `control-catallax.py run` with `--source recipe`, `--source field` and `--source compare`, using distinct output directories.

- Recipe must answer and write `answer.step`.
- Field must answer and write its own `display.obj`. Check the actual file, face count and display bound.
- Compare must report `comparison.verdict: agree`. Its `import_agreement` must contain 125 matching B-rep/field containment samples, both volume enclosures and a volume difference upper bound no larger than the allowed bound. A missing agreement object or a typed refusal does not pass.

Repeat field and comparison requests with `--workers 1`, `2` and `4`, separate outputs and `--fresh`. Compare manifests and work through the driver's `compare` command; also byte-compare the field OBJ files. Replay each exported `question.json` and require the same outputs and original allowance.

`cargo test --release -p catallax-compare --test import_step` generates public box, cylinder and 45-degree cone source fixtures and checks names, queries, volume, comparison, nonconvex planar membership and refusals.

## Coverage

Source admission, repairs and placement are shared with the B-rep. The compiled field keeps finite analytic definitions; its queries use no B-rep model. Current domains are polygonal planar parts and one cylindrical or 45-degree conical wall between perpendicular circular disks. Other curved trims return `field.import-domain`. Spheres and tori retain `field.import-unsupported`; spline faces retain `field.import-spline` with M10. Cone source angles retain the shared importer's degree-unit restriction. These refusals are coverage gaps, not impossible-geometry claims.
