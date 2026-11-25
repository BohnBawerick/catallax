Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# B02: retain a closed seam range's source-side parameters

P2 correction at reviewed head `32e460b40e6918b0c42d1666446217cd7ffb67b6`:
the native owning UV identity is now guarded through source-side pretrim,
initial clipping and opposite-side reclip. If trim enters short-span processing,
identity is not established and the seam path returns a typed refusal instead
of advertising zero owning drop with a changed parameterization. The reported
nonuniform helper mismatch was executed and confirmed; see
[before/after evidence](../../../../evidence/phase1/B02/P2/VALIDATION.md).
The same repair uses fixed-array chunk iteration in the existing 128-case test,
with an empty-remainder assertion and no lint suppression.

A rational cylinder cut by its tangent plane along the closed seam now returns one tangent range over the requested clipped domains, in either axial direction. Previously the full-period case returned `NotPorted(Intersect)` although either half-period returned one component. The repair retains the source iso-boundary's UV identity in `surface_trace.rs`; it does not port general ambiguous seam inversion.

**Partial B02 evidence:** candidate results and typed refusals are deterministic at actual 1, 2 and 4 validation workers. Application work equality is unmeasured because this direct candidate API exposes no request work counter. This is not full B02 work qualification, an application acceptance proof, or production parallelism.

## Baseline and consumer

Frozen main: `3a2f4e0c65feae548f4980d5c520a0c330acf6fb`, observed using `git ls-remote origin refs/heads/main`, fetched and checked out on 2026-11-14 UTC. Branch: `phase1/B02-seam-start`. The [merged P04 map](../P04/PORT-CONSUMERS.md) identifies the surface-intersector side-range call to `DropAndTrimCurve` as a seam-refusal risk for Phase 1.3 and indirectly 1.5. The current call remains in [surface_trace.rs](../../../../crates/catallax-trace/src/surface_trace.rs). [Trace candidates](../../../../crates/catallax-intersect/src/pair/trace_candidates.rs) run `Intersector::do_intersection` before independent complete-component selection/proof. B01's projected-drop no-consumer disposition does not apply to this true-branch range caller.

No named application recipe was replayed. This establishes a public candidate-level witness on that production route, not attribution of any specific housing or loft request. P04/B01 branches and evidence, other jobs and the P02 hold are preserved. No applicable `AGENTS.md` or `.agents/skills/SKILL.md` was found in the checkout/shared guidance. The existing verification skill was read for context; CLI qualification was not invoked.

## Public witness and repair

[The regression](../../../../crates/catallax-trace/tests/b02_seam.rs) constructs the radius-one rational cylinder about z using the nine circle control points and alternating weights `1` and `sqrt(1/2)` already described by the public [drop-curve generator](../../../../crates/catallax-trace/tests/fixtures/drop-curve/generate.py). The cylinder has degree (2,1), four rational quarter-circle spans and axial parameter z in [0,2]. The plane is x=1 with UV=(y,z). Cylinder axial domain is [0.25,1.75]; plane domain restricts z further to [0.5,1.5]. The full periodic cylinder domain is [0,4].

The cylinder's u-min and u-max sides are the same physical ruling. The existing long-range branch trims that side, then asks the general drop to invert it onto both surfaces. Inverting onto its own closed surface produces the two seam parameter solutions and hits the documented refusal. The curve already carries an unambiguous source-side identity; the search discarded that information before calling the drop.

The private range helper now receives `(surface index, side index)`. Only when the target is that source surface and the side's fixed axis is closed does it build the known UV line: `(fixed side u, t)` or `(t, fixed side v)`, with the trimmed curve's own interval. It checks endpoints against the requested source domain. The other surface's drop, trimming, distance filter, classification and existing duplicate removal stay in place. For trims that do not enter short-span processing, the owning side's drop distance is zero by construction; it is not a newly accepted numerical projection. The later independent Catallax proof boundary is unchanged.

This bounded construction avoids choosing an arbitrary seam solution. It does not change `surface_drop_curve.rs`'s general two-solution refusal, projected-drop refusal or tracer fallback. It does not add a general closed-surface seam walker, seam jump continuation or donor-parity claim for new positive outputs. The established source iso-curve identity justifies this range construction without C++ oracle regeneration.

## Executed validation

The final regression has 128 cases: 2 closed-axis choices (u/v), 2 axial orientations, 2 surface orders, 4 domain controls, and 4 periodic shifts. Domain controls retain the full period, retain the min-seam half, retain the max-seam half, or exclude the seam. Shifts translate the periodic knot vector and requested periodic domain by 0, 4, -4 and 8; rational geometric controls remain the same. This is parameter-period translation, not a claim about every geometric rigid translation.

Positive cases require exactly one tangent component, physical endpoints at z=0.5 and z=1.5, the expected reversed endpoint order when the axial controls reverse, and both UV control polygons inside their requested domains without an allowance. Seventeen samples check each curve and surface image with a test-only residual threshold of 1e-12. This threshold is not a changed production tolerance or proof certificate. Floating-point spline evaluation of an exactly constant boundary coordinate can round slightly outside it; the stored UV control coordinates are still checked with zero allowance. Seam-excluding controls require an empty result.

Period-translated cases compare the complete stored output bits after subtracting only the known period shift from the owner's UV control coordinate. They compare degree, dimension, knot/control counts, every knot and homogeneous control coordinate, component order, kind, deviation and UV presence. Worker comparisons additionally compare the **unnormalized** output bits and exact typed refusal strings for each case. Every case also checks that arbitrary general seam drops with both values of `on_surface` retain `NotPorted(Intersect)`.

Actual scoped OS threads synchronize at a barrier. Case index i belongs to worker `i % workers`; results are reassembled in original index order. Partitions are [128], [64,64], and [32,32,32,32]. All partitions are nonempty. Each case owns its intersector and surfaces; there is no shared candidate mutation. This parallelizes validation cases only. No production scheduling, request allowance, work charge, caps, tolerance or saved expectation changed.

| Check | Observed result |
| --- | --- |
| New translated/clipped regression on frozen baseline source | Exit 101; first full-period case returns `NotPorted(Intersect)` |
| Same regression after restoring candidate source | Exit 0; all 3 test functions pass; full matrix and worker comparisons pass |
| Existing `curve_surface` integration suite | 6 passed |
| Existing `drop_curve` integration suite | 3 passed, including the two saved refusal rows |
| Existing `intersector` integration suite | 6 passed, including requested-domain, saddle and singular-inside-branch controls and the saved differential gate |
| Static checks | Rustfmt on the two changed Rust files; relative links, public machine-path audit, source hashes and `git diff --check` |

The finalized regression was also run against the frozen baseline implementation by temporarily restoring that source file, then restoring the candidate in a `finally` block. The sanitized [baseline log](../../../../evidence/phase1/B02/baseline-regression.txt) preserves the observed failure; its original is retained privately. Earlier exploratory six-case output showed the full seam refusal, both half-domain successes and an already-successful transverse ring cut. Earlier assertion development encountered ordinary floating-point endpoint/evaluation roundoff; final tests distinguish stored-domain bounds from sampled image residuals. No production tolerance was changed to satisfy those checks.

[Validation record](../../../../evidence/phase1/B02/VALIDATION.md) records portable commands and actual summaries; [artifact hashes](../../../../evidence/phase1/B02/artifact-digests.json) identify the final public files. Test outcomes are Rust observations, not fresh donor-oracle evidence.

## Remaining requirements and cleanup

Application-level work, the one-request allowance outcome, Boolean/STEP success and independent import were not measured. The direct candidate layer has no application work counter; do not equate output size, sample count or test time with deterministic application work. No new counter, allowance or reset was introduced. General free seam starts still refuse where appropriate. The repair resolves one bounded seam-range family; it does not qualify every closed surface or implement both directions of the general tracer's seam continuation. Both axial directions here are validated through reversed surface controls and opposite component endpoint order.

No C++ execution, private/raw corpus, companion repository, sanitizer/memory-fault probe, hosted gate read/retry, CI launch, Result/build adapter or PRD change occurred. No memory-fault work was needed. No PR was created, merged or deployed. Disposable source-restore scratch was removed after preserving the baseline log; reusable build cache, checkout, user work and required evidence are retained. Exact pushed head and remote verification are provided in the handoff.
