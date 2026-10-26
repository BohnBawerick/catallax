## 2026-10-26: unproved B-rep preview mesh

Pieces: display mesh, application entry

Implements the early mesh and selection delivery lane in [ARCHITECTURE.md, "Preview and deferred proof: decision of 2026-10-26"](../ARCHITECTURE.md#preview-and-deferred-proof-decision-of-2026-10-26), approved in https://github.com/BohnBawerick/catallax/pull/129.

### API

`catallax_mesh::preview_mesh_body(model, body, &PreviewOptions, decider)` returns `Result<PreviewMesh, Refusal>`, never `Answer<Mesh>`. `PreviewMesh` contains binary64 positions, face-labelled triangles, `PreviewFace` picking ranges and `PreviewEdge` polylines. Face and edge IDs refer to the retained input model and its origins. It contains no certificate, marks presented as proof, deviation bounds or quality bounds. It has no export writer or conversion to certified `Mesh`.

`PreviewOptions` takes the same chord and angular targets as certified display, plus `max_vertices`. `PreviewOptions::from(&MeshOptions)` copies only those three controls. Chord and angle are visual sampling targets, not guaranteed limits. Certified edge-length and aspect-ratio controls do not belong to this API. Resource or precision failures return no partial preview.

`catallax_app::preview_body(&Built, &PreviewOptions, decider)` returns `PreviewDisplay { mesh, names }` from an existing retained build. The name table contains complete face origins and result edge/vertex names for that same revision. `catallax_app::build_preview(recipe, options, decider)` returns `PreviewApplicationBody { built, preview }` without requesting a certified display mesh. SolidDraft-format callers can first use `soliddraft::convert`; incremental callers can pass `checkpoint.built()` to `preview_body`.

Charge `preview_body` to the preview/display allowance. It uses an active caller allowance unchanged. Without one, the existing `DEFAULT_BUILD_WORK` limit applies. `build_preview` shares one allowance across construction and preview. Neither entry starts background work, resets the decider deadline, spawns numeric threads, duplicates an allowance or reserves proof work. A later proof request uses the retained `built.model` and `built.result`; the preview does not authorize it or enlarge its allowance.

### Sampling and unchanged guarantees

The preview samples retained edge traces and face definitions. Analytic ruled strips share sampled angular columns across their end trims. Other analytic trims use the existing periodic chart triangulator with sampled refinement. The existing retained spline-prism family uses sampled paired rings and cap fans. No preview path calls `mesh_body` or certifies facets. Preview topology, watertightness and approximation quality are unproved.

`mesh_body`, `Answer<Mesh>`, every certified bound, `Mesh::write_obj` and its `proven_chord_error` labels, certified watertightness checks, STEP and all other export guarantees remain unchanged. There is no automatic export fallback to preview data.

### Development measurements

Run the development-only benchmark from the repository root, with the JSON output outside the checkout:

```sh
taskset -c 2 cargo bench -p catallax-app --bench preview_mesh -- /tmp/preview-mesh.json 3
```

It builds each public fixture once, warms both mesh paths, alternates their order, and uses fresh deciders with identical 20,000,000-work allowances and 10-second deadlines for each display request. Rows separate elapsed time, deterministic work, output sizes and refusals. Fixtures include the public fillet and countersink reproductions, the [SolidDraft chamfer scenario](https://github.com/BohnBawerick/soliddraft/pull/30) through its distance edit with the pocket retained, P2, P1 and P3. The optional third argument selects one fixture by name.

Eight diagnostic points per preview triangle compare sampled positions and normals to retained analytic carriers. They are explicitly labelled `diagnostic_not_a_proof`. These samples do not prove trim coverage, Hausdorff distance, watertightness, deviation bounds or export eligibility. Timings include application picking names on the preview side and certified meshing alone on the certified side. Building the body, sampled diagnostics, JSON and any later proof are outside those display timings. Shared-machine timing remains diagnostic; compare warmed, pinned pairs rather than a single sample.

Revision lifecycle states, stale-preview replacement, deferred body construction proofs, field preview and certificate reuse remain follow-ups. A caller still owns revision identity and must not resolve numeric picking handles against another model.
