# Unproved GPU field preview

This crate holds the field's unproved graphics-card routes: a batch point evaluator, a coarse first paint of the field display, and the throwaway F9 viewer. Nothing here is proved. Values, meshes, pixels and timings never feed an export, a proof or a proven measurement.

## Batch evaluator

`batch::Evaluator` evaluates the field's unproved WGSL program in wgpu compute shaders. It takes arbitrary points or a regular grid and returns `Samples`: one value per point in input order, negative inside and positive outside, with the side, a distance estimate and the plan's true-distance flag. When that flag is false the magnitude is a Boolean residual, not a distance. Every batch records the hardware that made it.

The processor preview, `Plan::unproved_preview_batch`, is the fallback. It runs when no device is given, when no adapter or device exists, or when the plan has no WGSL translation; the reason is recorded in `Hardware::Processor`. A shader that fails validation refuses with `preview.shader` instead of falling back.

Graphics-card arithmetic is binary32. Arbitrary points are rounded to binary32 before evaluation. Grid points are cell centres, `(origin + index + 1/2) * 2^e` with indices below 2^23, so they are exact in binary32 and binary64 and both routes evaluate identical points. The half step keeps samples off faces at whole multiples of the pitch, where a zero value would leave the side to the last bit. Each invocation writes one value, with no reduction or atomic, so repeated runs and every dispatch size give the same bits on one device. Processor values are the same bits at any worker count. A non-finite graphics-card value refuses the whole batch with `preview.gpu-non-finite`; no partial batch is kept.

## First paint

`first_paint::first_paint` samples a plan on a power-of-two grid covering its exact material bounds plus two pitches, then extracts a surface-nets mesh on the processor. A grid above the caller's point limit refuses with `preview.grid-size` before any evaluation; `MAX_POINTS`, 2^28 points, keeps the binary64 values near 2 GB. A point is inside when its value is negative. The `SampledMesh` has outward-wound triangles, a count of open edges at the grid border, and its enclosed volume. It has no face ownership, writer, deviation bound or conversion to a certified mesh. The grid pitch is a visual setting, not an achieved error bound.

`catallax-first-paint` checks the first paint against the processor preview on P1 to P4:

```sh
cargo run --release -p catallax-preview --bin catallax-first-paint -- "$ARTIFACTS"
cargo run --release -p catallax-preview --bin catallax-first-paint -- "$ARTIFACTS" --cpu --workers 4
```

The output directory must be outside the repository. By default it samples pitches of 0.5, 0.25 and 0.125 mm (`--pitch -1 --pitch -2 --pitch -3`); `--part P3` limits the parts and `--max-points N` changes the grid limit. For each part and pitch, `report.json` records:

- graphics-card and processor sample times, a repeated graphics-card run and a run split into dispatches of 100,003 points, each compared bit for bit;
- processor values at 1, 2 and 4 workers on the coarsest pitch, compared bit for bit;
- side agreement with the processor preview at every grid point, the largest value difference, and the first 64 mismatching points with both values;
- the graphics-card mesh and the processor-sampled mesh: vertices, triangles, open edges, volume, whether the triangle lists are identical and the largest vertex shift;
- distances from the graphics-card mesh's vertices to the processor chart preview (`catallax_field::display::preview_plan`, under the 20M application allowance) and back, with that preview's cells, unresolved cells, work and any refusal.

It also evaluates 1,000,000 seeded points through the arbitrary-point route and compares them at their binary32 coordinates. A run where graphics-card values change between runs or dispatch sizes, or processor values change with the worker count, exits with `preview.nondeterministic` after writing the report. Mismatches and timings are information only.

## Throwaway F9 viewer

This throwaway F9 viewer renders the public P3 recipe at 960 by 720, then rebuilds it after changing the junction radius from 2 mm to 3 mm. It writes both frames, per-pixel face IDs, hit points and reports outside the repository. It has no export or proven measurement output. SolidDraft first-display integration is the remaining half of M15.

```sh
cargo run --release -p catallax-preview -- "$ARTIFACTS" --workers 4
cargo run --release -p catallax-preview -- "$CPU_ARTIFACTS" --cpu --workers 4
```

Set the output variables to directories outside the repository. Without `--workers`, processor work uses the available cores. Each GPU invocation owns one pixel. A missing adapter or unavailable device selects the processor preview; `--cpu` exercises that route deliberately. A shader failure, unfinished ray or missing reference returns a typed refusal.

The WGSL field translator reads the compiled tape and its cached constants. It translates stock solids, finite and through extrusions, full and partial revolutions, analytic cells with supported surfaces, and admitted planar, polyhedral, round, annular and oblique draft cells. Other leaves and constants outside binary32 range refuse the translation. Per-point arithmetic is ordinary binary32 on the GPU. A 20,000-point check for each frame compares GPU values with the processor preview before rendering. These checks give sampled evidence, not a field error certificate.

The ray marcher advances by the field residual and tests the sign transition. A Boolean residual can be zero on a buried support, so zero alone never becomes a visible hit. Each ray has 16,384 steps and each crossing has 16 bisections. Unfinished rays refuse the frame.

The GPU also computes an unproved face ID for every hit. This ownership evaluator is a specialization for P3's thirteen closed analytic patches: planes with finite trims, the boss cylinder, its disk, the junction torus quarter and the four trimmed rim quarters. It minimizes distance to the finite patches rather than picking a winning Boolean instruction or an untrimmed support. Rim minimization uses 32 angular intervals and 28 refinements. Recipe names come from the compiled field registry. The processor fallback evaluates the same patches with ordinary binary64 arithmetic. This specialization is not a replacement for F8's general certified query.

Both frames and their ID buffers are saved before exact selection checks. Each frame then receives 1,024 unique seeded samples drawn from all surface pixels, with no ownership-border exclusion. The fixed seed is `0xf920261111`. Each sampled binary32 hit is converted to its exact rational coordinates and queried through F8's `Prepared::nearest_face`, with unchanged default controls and the existing F8 regression allowance of 100,000,000 work per query. The report gives agreement over all samples, proved and refused counts, and every non-agreeing sample's pixel, point, preview ID/name and exact F8 answer. Refusals remain explicit. A proved mismatch fails the run after both frames' checks finish.

Each pixel has a little-endian u32 entry: zero is background, and other IDs index the accompanying recipe-name table. IDs and hit positions remain unproved. They never become a proven measurement or an export input.

The independent reference builds and checks the B-rep from the same recipe. It rasterizes its default 0.01 mm chord-tolerance mesh with the same orthographic camera. The report records the symmetric Euclidean Hausdorff distance between the two pixel outlines. A distance above one pixel fails. This is an outline check for this camera, not a global distance certificate or an export check.

For each revision the output contains:

- `<name>.ppm`: face-colored preview with a visible unproved-preview label.
- `<name>-brep.ppm`: B-rep reference mask.
- `<name>.face-ids.u32le`: 960 by 720 row-major face IDs.
- `<name>.hits.f32le`: row-major four-component binary32 records, with xyz hit coordinates and status zero for background or one for a hit.
- `<name>.names.json`: IDs and recipe names.
- `<name>.f8-samples.json`: all exact sample answers, query work and every mismatch.
- `report.json`: backend, dimensions, radius, outline error, sampled evaluator differences, ownership agreement and elapsed times. Pending sample checks are labelled as pending.

Frame, ownership and exact-reference times are information only. Raw evidence is retained privately.


## General named-face picking

`picking::Picker::new(&chart_preview)` snapshots the field chart preview's named
triangles. `picker.points(&points, gpu, workers, max_tests)` returns face handles,
triangle-distance estimates, backend details, missing-cell counts and omitted
collapsed-triangle counts. It accepts arbitrary parts and poses. It needs no P3
recipe constants or B-rep build. The caller resolves handles against the same
revision's field registry.

These selections are unproved. A nearest triangle in a picture with gaps need
not be the nearest defined face. Equal minima select the first triangle in the
picture. The processor uses binary64 after rounding input coordinates and
triangles to binary32, and returns the same selections and distances at every
worker count. The GPU uses binary32. `max_tests` limits the point-triangle
product before work starts; typed refusals include the counts. The
[interface record](../../docs/interface-changes/2026-11-11-unproved-gpu-preview.md#general-named-preview-picking)
records the other limits and codes.

The surface-nets first paint remains unnamed. Use its named chart preview when
available, or `Prepared::nearest_face` for a certified pick with explicit
ambiguity and distance bounds. The F9 binary retains its P3 specialization and
its existing pixel checks. General picking is a reusable library route for the
application, whose first-display integration is separate.

Run the qualification checks with:

```sh
cargo test --release -p catallax-preview --test shader_coverage -- --nocapture
cargo test --release -p catallax-preview --test picking -- --nocapture
```

Set `CATALLAX_REQUIRE_GPU=1` to fail a qualification run when no adapter is
available. The shader suite also validates every program without a device and
without proof work. It covers P1-P4, drafted P1, the admitted P5 radius-1/2 edit,
draft controls, partial turns and a through cut. Authored P5 at radius 2 retains
`fillet.unsupported-overflow` before translation; the edit is a separate case.
