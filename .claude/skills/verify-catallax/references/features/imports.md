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

`cargo test --release -p catallax-compare --test import_step` generates public box, cylinder, cone (45 degrees, other radian and degree angles, and pointed), sphere (including the `VERTEX_LOOP` form) and torus source fixtures. It checks names, queries, volume, comparison, recorded cone tolerances, nonconvex planar membership and refusals.

## Coverage

Source admission, repairs and placement are shared with the B-rep. The compiled field keeps finite analytic definitions; its queries use no B-rep model. Current domains are polygonal planar parts, one cylindrical or conical wall between perpendicular circular disks, a pointed cone with one disk, and a complete sphere or torus. Cones at any semi-angle are built through their rims; the recipe's import record lists each such cone in `face_tolerances`, at most 0.01 mm. Other curved trims, including trimmed spheres and tori, return `field.import-domain`. Spline faces retain `field.import-spline` with M10. A torus whose axis has irrational length returns `field.import-domain`. These refusals are coverage gaps, not impossible-geometry claims.
