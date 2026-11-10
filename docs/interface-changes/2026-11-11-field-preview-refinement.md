# Progressive field preview

`Revision::construct_progressive` adds caller-scheduled refinement to the native
field construction route. Existing `construct`, certified display, measurements
and export calls keep their contracts.

```text
construct_progressive(recipe, options, proof_options, allocation)
  publish 0.1 mm first paint, unproved
  return RefinementJob
background: RefinementJob::refine
  publish 0.01 mm refinement, unproved
  return ProofJob
background: ProofJob::prove().publish()
  confirm the same revision or attach its typed refusal
```

The application chooses when to schedule each background job. Catallax starts no
thread. Both paints retain the same immutable recipe. A newer edit invalidates
the old job, so a late refinement cannot replace the newer display. Dropping a
job postpones its remaining work and leaves the displayed mesh unproved.

`ProgressiveOptions::first_paint_work` splits `WorkAllocation::preview` between
the paints. `limits` supplies each paint's existing cell and vertex limits. The
method overrides its visual tolerance with 0.1 mm and 0.01 mm respectively.
Refinement includes field recompilation. All three work caps are reserved from
the caller's allowance without increasing the total. Unused work does not move
between phases. `WorkReport::preview` accumulates actual work from both paints.
An incomplete display, exhausted allowance or unsupported definition retains
its typed refusal. Failed refinement marks the affected revision failed.

Both paints use the existing field mesher and are `Status::Unproved`. Their
tolerances are visual sampling targets, not proven distance bounds.
`PreviewMesh` has no certified mesh
conversion or export writer. Export and proven measurement consumers still
require the confirmed exact body and their own scoped proof. The coarse first
paint is a 0.1 mm display target. The refinement uses today's 0.01 mm mesh.

For diagnostic runs, build `cargo build --release -p catallax-field --bench
preview -j 2` and run the resulting `preview` executable with `RECIPE OUTPUT_DIR
certified|preview TOLERANCE_DENOMINATOR`. A denominator of `10` requests the first
paint's target, and `100` the refinement target. The executable fixes the work
allowance at the verification contract's 100M units. Optional
`CATALLAX_PREVIEW_DEADLINE_SECONDS` limits a diagnostic run and retains its typed
deadline refusal. Timings include field compilation and display. They are
informational shared-machine numbers, with no unloaded-machine guarantee or
M26 speed qualification.

Successful runs retain `mesh.bin`, a little-endian binary64 vertex and triangle
payload, and `report.json`. Certified runs also retain the actual `display.obj`.
The payload starts with a u64 vertex count, then three binary64 coordinates per
vertex, followed by four u64 values per triangle: three indices and a face tag.
Compare payload and OBJ bytes across revisions, and inspect unresolved-cell and
refusal fields before treating a run as complete. Raw evidence retained privately.
