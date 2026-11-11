# Unproved GPU field preview

This throwaway F9 viewer renders the public P3 recipe at 960 by 720, then rebuilds it after changing the junction radius from 2 mm to 3 mm. It writes review frames, hit points, face IDs and a report outside the repository. It has no export or proven measurement output. SolidDraft first-display integration is the remaining half of M15.

```sh
cargo run --release -p catallax-preview -- "$ARTIFACTS" --workers 4
cargo run --release -p catallax-preview -- "$CPU_ARTIFACTS" --cpu --workers 4
```

Set the output variables to directories outside the repository. Without `--workers`, processor work uses the available cores. GPU invocations each own one pixel. A missing adapter or unavailable device selects the processor preview; `--cpu` exercises that route deliberately. A shader failure, unfinished ray, exhausted ownership allowance or missing B-rep reference returns a typed refusal instead of a passing report.

The WGSL translator reads the compiled field tape and its cached constants. It translates stock solids, finite extrusions, full revolutions and analytic cells with supported surfaces. Other leaves and non-finite binary32 constants refuse the translation. Per-point arithmetic is ordinary binary32 on the GPU. A 20,000-point check for each frame compares GPU values with the processor preview before rendering. These checks give sampled evidence, not a field error certificate.

The ray marcher advances by the field residual and tests the sign transition. This matters because a Boolean residual can be zero on a buried support. A zero residual alone never becomes a visible hit. Each ray has 16,384 steps and each crossing has 16 bisections. Unfinished rays refuse the frame.

F8 ownership runs on the processor after GPU readback. It queries the exact rational representation of each binary32 hit point. A proved open ownership ball can supply IDs to nearby points inside half its radius. Fixed 16-row tiles keep the cache, query order and logical work independent of worker count. Other points call `Prepared::nearest_face` with its unchanged default controls and the existing F8 regression allowance of 100,000,000 work per query. The viewer never infers an owner from a winning Boolean instruction. Every pixel has a little-endian u32 entry: zero is background, `4294967295` is an F8 ambiguity refusal, and other IDs index the accompanying recipe-name table. IDs name faces but do not prove the GPU hit position. Refused pixels are red.

The independent reference builds and checks the B-rep from the same recipe. It rasterizes its default 0.01 mm chord-tolerance mesh with the same orthographic camera. The report records the symmetric Euclidean Hausdorff distance between the two pixel outlines. A distance above one pixel fails. This is an outline check for this camera, not a global distance certificate or an export check.

For each revision the output contains:

- `<name>.ppm`: face-colored unproved preview, with a visible unproved-preview label and a PPM header label.
- `<name>-brep.ppm`: B-rep reference mask.
- `<name>.face-ids.u32le`: 960 by 720 row-major face IDs.
- `<name>.hits.f32le`: row-major four-component binary32 records, with xyz hit coordinates and status zero for background or one for a hit.
- `<name>.names.json`: IDs and recipe names, including the refusal sentinel.
- `report.json`: backend, dimensions, radius, outline error, sampled evaluator differences, query/reuse/refusal counts, query work and elapsed times.

Frame and ownership times are information only. Raw evidence is retained privately.
