## 2026-10-26: unproved field preview

Pieces: field display

Implements the field half of the preview lane in [ARCHITECTURE.md, "Preview and deferred proof: decision of 2026-10-26"](../ARCHITECTURE.md#preview-and-deferred-proof-decision-of-2026-10-26), following the B-rep pattern in [2026-10-26-unproved-preview-mesh.md](2026-10-26-unproved-preview-mesh.md).

### API

`catallax_field::display::preview_plan(plan, &PreviewOptions, decider)` returns `Result<PreviewMesh, Refusal>`. `PreviewMesh` holds binary64 positions, triangles tagged with field face handles (origins stay resolvable through the plan), the visited `cells` and `unresolved_cells`. It carries no `max_deviation`, certificate or quality bound, and has no OBJ writer or conversion to the certified `Mesh`. `PreviewOptions { tolerance, max_vertices, max_cells }` takes the certified controls; `PreviewOptions::from(&Options)` copies them. The tolerance is a visual sampling target, not a guaranteed limit.

The preview reads the same operation definitions and uses the same field charts, trim predicates and cell classification as `mesh_plan`. It differs only where the certified request would refuse a cell: a cell that exhausts the cell limit, stays unresolved or fails the rounding check is counted in `unresolved_cells` and left empty. Nothing fills a gap, and no B-rep geometry, sample or hull is used. A nonzero `unresolved_cells` means the preview is incomplete, and a zero count is not a certificate either.

The caller's allowance is charged unchanged. An exhausted operation allowance or deadline, the vertex limit, an unsupported chart and lattice leaves still return a typed refusal with no partial preview.

### Unchanged

`mesh_plan`, `mesh_plan_reported`, `Options`, `Mesh`, `Mesh::write_obj`, every limit and every refusal code of the certified request, and its bytes. The `field` verification adapter is not changed and still waits for the certified mesh; wiring the preview into it is a follow-up once the adapter's open change lands.

Development timing is in the ignored test `preview_against_certified_time_on_public_fixtures`. Where the certified request succeeds the preview does the same work, so it is not faster there; the preview exists for requests the certified display refuses. Shared-cloud timing is diagnostic only.
