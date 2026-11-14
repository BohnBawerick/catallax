# B02 P2: native parameter identity through seam-range trims

Repair base and verified starting remote/PR head:
`32e460b40e6918b0c42d1666446217cd7ffb67b6`, branch
`phase1/B02-seam-start`, open [PR #607](https://github.com/BohnBawerick/catallax/pull/607).
The independent [P2 review](https://github.com/BohnBawerick/catallax/pull/607#pullrequestreview-5907767750)
was read before editing. The clean local worktree, remote branch and PR head
matched. No rebase or merge from current main was performed; N02 and other work
remain separate. Applicable checkout guidance was checked; no AGENTS.md or
.agents/skills/SKILL.md was present. The verification skill was read for context,
not invoked as an application/STEP qualification workflow.

## Executed before/after helper witness

The test-only private-helper module is
`crates/catallax-trace/tests/support/b02_parameter.rs`. It constructs a closed-u
radius-one rational cylinder with degree-one axial knots [0,0,.5,1,1] and heights
[0,1,1+1e-7]. A supplied source-side range [.125,1] is passed to the private
`side_range_curve` helper against a tangent plane. These are public constructed
inputs, not donor calls or a claim that the public solver emits this range.

Before the production edit, the four new helper tests executed against the
reviewed implementation: one uniform-control test passed and three tests failed
on interior correspondence. The source prediction is confirmed at the helper:

```text
nonuniform=true pretrim=false reclip=false swap=false t=.56
spatial_z=0.6228571925714287 lifted_z=1.000000012 max_gap=0.42187495781249984
nonuniform=true pretrim=true reclip=false swap=false t=.56
spatial_z=0.6228571925714287 lifted_z=1.000000012 max_gap=0.42187495781249984
nonuniform=true pretrim=false reclip=true swap=false t=.56
spatial_z=0.6228571925714287 lifted_z=1.000000012 max_gap=0.42187495781249984
```

The initial range is [.125,1]. The source-pretrim control supplies source domain
v=[.125,1]. The opposite-drop-reclip control starts from the full [0,1] range
and clips the plane's z domain below at .25; the opposite drop restricts the
varying interval before the spatial reclip. The test checks 33 interior/end
sample locations through the owning UV, in addition to printing t=.56.
No test asserts that these samples establish a whole-curve geometric proof.

After the repair, all four helper tests pass. The nonuniform controls return
`SolverRefusal::NotPorted(Intersect)` with no emitted component at source pretrim,
initial clip, and opposite-side reclip. Reclip is exercised with either surface
order. The bounded uniform control uses the same knots and heights [0,.5,1];
all six combinations of surface order and clipping stage still produce exactly
one component with max sampled owning gap 0 and z=.56 at t=.56.

This is a conservative refusal, not an inferred parameter remap. Whenever a trim
needed for the native owning seam identity reports `short_spans`, it may have
removed segments and redistributed knots, so the port refuses even if a
particular removal ultimately changed nothing. Typed underlying trim failures
also propagate on this seam path. Legacy non-seam trimming is unchanged.
The guard covers source-side construction, initial spatial clipping and both
spatial/UV back-trims after the opposite-side drop. No tolerance/cap or general
curve-drop implementation is changed.

## Commands, controls and limits

Portable command forms from the repository root:

```sh
cargo test --offline -p catallax-trace --lib b02_parameter -- --nocapture
cargo test --offline -p catallax-trace --lib --test b02_seam --test intersector --test drop_curve
ci/run clippy -p catallax-trace
rustfmt --edition 2024 --check crates/catallax-trace/src/surface_trace.rs crates/catallax-trace/tests/support/b02_parameter.rs
git diff --check
```

Actual commands selected existing authorized Cargo/Rustup roots through
local environment variables; there was no installation or permission change.
Before/after observations above are transcribed from actual tool output, not
purported raw stdout files. All four helper tests passed after repair.

The same authorized follow-up replaces `cases.chunks_exact(4)` with
`let (families, remainder) = cases.as_chunks::<4>();`, asserts an empty remainder,
and iterates the fixed-array families. Local official all-target Clippy first
reproduced `clippy::chunks_exact_to_as_chunks` under warning denial, then passed
after the correction. No lint suppression was added. The installed compiler
reports `rustc 1.100.0 (0ef1bc386 2026-11-09)`; the fixed-array API is supported.
No hosted job was read or rerun. The 128 inputs and complete comparisons are
unchanged; the final matrix run is recorded below.

The original B02 captures remain unchanged and their original metadata is
preserved in the old commit and authorized private scratch. The root B02 digest
manifest identifies current source/docs and the retained historical capture.
The prior B02 128-case matrix checks complete candidate bits at actual 1/2/4
validation workers, separately from these serial helper witnesses; it does not
measure application work. No new production parallelism, application work
qualification, public end-to-end occurrence, Boolean/STEP acceptance, donor
parity or general seam continuation is claimed. Independent re-review is left to
the supervisor; no comment, new PR, merge, hosted CI, B03 work, C++ or private
input is part of this repair.

## Final focused results

- Combined focused command: exit 0, 9 library tests (including the four helper
  witnesses), 3 B02 seam tests, 3 drop-curve tests and 6 intersector tests pass.
  Saved expectations are unchanged. The intersector gate completes in 75.61 s;
  this is a shared-machine observation, not a benchmark.
- Final post-iterator B02 run: exit 0, 3 tests pass; 128 cases at actual worker
  counts 1, 2 and 4 retain identical complete outcomes and refusal records.
  Partitions remain [128], [64,64], [32,32,32,32]. Application work is unmetered.
- `ci/run clippy -p catallax-trace`: exit 0; all targets, warnings denied.
- Rustfmt checks all three changed Rust files; diff hygiene passes.
- Original baseline capture bytes are unchanged; public path/link and digest
  audits pass. No N02 or unrelated branch files are modified.
