## 2026-11-11: field display stage clocks

Pieces: catallax-field, verification CLI, verification driver
Change: `catallax_field::display::Cost` adds `elapsed: Duration`, the observed
wall time of the operations its counters describe. `mesh_plan_reported` fills
it for `total`, `common`, every support stage and every face query; `common` is
the remainder of `total` after the support totals. Field verification
`display_diagnostics` adds `wall_ms` to each cost record and a `compile` record
(`decisions`, `work`, `wall_ms`) for the request's field compilation. The
driver's `comparable` drops every `wall_ms` inside `display_diagnostics`, and
its new `field-stages` command prints a `parts` report's stage table, clocks
beside work.
Why: PLAN M13 F2. A display wait can only be shortened once its stages are
measured where the work counters already are.
Migration: code that builds `Cost` with a struct literal adds `elapsed` or uses
`..Default::default()`. Code comparing `Cost` or `Diagnostics` values from two
runs compares the counters, not the whole value. Callers comparing field
manifests exclude `display_diagnostics` clocks as they exclude
`first_mesh_ms`. Display bytes, work, decisions, cell counts and outcomes are
unchanged; field manifests gain the additive `compile` record.
