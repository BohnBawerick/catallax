# Unproved GPU field preview

This throwaway F9 viewer renders the public P3 recipe at 960 by 720, then rebuilds it after changing the junction radius from 2 mm to 3 mm. It writes both frames, per-pixel face IDs, hit points and reports outside the repository. It has no export or proven measurement output. SolidDraft first-display integration is the remaining half of M15.

```sh
cargo run --release -p catallax-preview -- "$ARTIFACTS" --workers 4
cargo run --release -p catallax-preview -- "$CPU_ARTIFACTS" --cpu --workers 4
```

Set the output variables to directories outside the repository. Without `--workers`, processor work uses the available cores. Each GPU invocation owns one pixel. A missing adapter or unavailable device selects the processor preview; `--cpu` exercises that route deliberately. A shader failure, unfinished ray or missing reference returns a typed refusal.

The WGSL field translator reads the compiled tape and its cached constants. It translates stock solids, finite extrusions, full revolutions and analytic cells with supported surfaces. Other leaves and constants outside binary32 range refuse the translation. Per-point arithmetic is ordinary binary32 on the GPU. A 20,000-point check for each frame compares GPU values with the processor preview before rendering. These checks give sampled evidence, not a field error certificate.

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
