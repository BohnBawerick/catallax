# R11 first step: refusal census, source gap and historical control

**11 November 2026: today's census is blocked by a missing current log.**
The requested interval is 2026-11-11 00:00–24:00 UTC. No attributable operation-log
bundle for that interval was supplied or found among the tracked project logs.
The current-day denominator and counts are **unknown**, not zero. This note
fulfils only the source audit and a historical control; it does not complete
[PLAN.md R11](PLAN.md#126-robustness-and-many-core-steps) or justify step two's
field-based surviving-face prototype.

Audit baseline: `25b7c2d1341db65e0b8ee494a0f4309d95259008`. No new failure
campaign was run. Only public committed material was used; no private inputs,
saved-answer corpus, companion repository or C++ oracle was accessed for counts.

## Eligible source and counting rule

The kernel log is the verifier's `catallax-operation-log/1`, embedded in
`operation_log` and also written as `operations.json`; these are two copies of
one log, not two executions. See the [log contract](VERIFY-INTERFACE.md#operation-log)
and [writer](../crates/catallax-cli/src/verification.rs). The schema has no UTC
execution timestamp or build identity at the log root. A dated run manifest,
binary/build identity and public input provenance must accompany it.

The available historical control is the public [22 October parts snapshot](checkpoints/2026-10-22-parts.json),
with provenance in its [checkpoint note](checkpoints/2026-10-22-verification.md):

- Documentary date: **2026-10-22**. UTC start/end timestamps are not recorded;
  do not infer a precise execution interval or recency from this date.
- Recorded build revision: `d98ac6a00552760399e2cd5f58b55a4a90179c0b`;
  binary SHA-256 `6dc829c48dc39381dd4557afb8b4ef78866142ae189c3f57894a9f91393b26ff`.
  These are historical provenance claims, not a fresh build verification.
- Snapshot SHA-256: `b3d1eef01215a9a679be1790badf3e3cb4ede49c6911fe63959cb5935fce7eee`.
- Unchanged recorded request allowance: 1,000,000 work units.
- Scope: 13 recipe records × 3 adapters (`recipe`, `field`, `compare`) = **39 logs**,
  containing **240 executed spans**. Only **2 spans**, in **1 unique recipe**,
  have operation exactly `shell`, `offset` or `thicken`.

Count only the immediate `cases[i][adapter].operation_log.operations` array.
Deduplicate by (snapshot digest, recipe SHA-256, adapter, kernel, sequence,
operation ID, operation name); identical duplicate keys must have identical
records or become an audit error. Root refusals and manifest copies add no events.
Distinct adapters are separate executions. Compare preflight outcomes and later
measurement/export failures are not attributed to unexecuted construction steps.
An empty execution prefix supplies no successful shell/offset/thicken evidence.

## Historical observed events — not today's failures

| Operation | Executed | Answered | Refused: unsupported scope | Proof/work budget | Invalid input/selection | Internal error/interrupted | Unclassified refusal |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Shell | 2 | 0 | 2 | 0 | 0 | 0 | 0 |
| Offset | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Thicken | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Total | 2 | 0 | 2 | 0 | 0 | 0 | 0 |

Both shell spans refuse before numeric decisions or work (zero in each).
They are observed conservative refusals, **not observed correctness defects**.
The unsupported-scope share is 2/2 selected executions, or 1/1 unique selected
recipe; no rate exists for offset or thicken (0/0). This tiny, obsolete sample
cannot estimate current capability or the frequency of surviving-face problems.

| Cause | Count | Code and representative public record |
| --- | ---: | --- |
| B-rep shell operation not yet implemented in the recorded build | 1 | `features.unsupported-operation`; snapshot JSON pointer `/cases/11/recipe/operation_log/operations/1`, kernel `brep`, ID `tray`, sequence 1 |
| Field shell operation not yet implemented in the recorded build | 1 | `field.unsupported-operation`; snapshot JSON pointer `/cases/11/field/operation_log/operations/1`, kernel `field`, ID `tray`, sequence 1 |

Both identify `P4-shelled-tray.json`, recorded recipe SHA-256
`94eda6ca26261257497ee8c2747f956e4b33806e100101111ece65d857d510fd`.
The snapshot records, not today's edited fixture, define the historical case.
The P4 compare adapter has an empty operation array: it contributes no third
shell event. Other 238 spans are outside the selected operation names.

## Current-source cause map — counts not measured

The following are classification aids at the audit baseline, not logged events
or fresh test results. Code names alone are insufficient to determine cause.

| Category | Source-supported distinction | Representative source/public control |
| --- | --- | --- |
| Unsupported construction or version scope | `distance.unsupported-envelope` can mean a missing normal-envelope route or legacy sharp-join semantics, not a numeric exhaustion | [sharp selection control](../crates/catallax-features/tests/sharp_thicken.rs), `version_three_keeps_the_earlier_sharp_selection_refusal`; [spline shell controls](../crates/catallax-features/tests/spline_cap_shell.rs), `a_spline_pocket_is_not_mistaken_for_a_complete_prism` |
| Invalid thickness or selected reference | `distance.non-positive-thickness` and `reference.lost-face` distinguish invalid/non-surviving selections from construction coverage | [dispatch](../crates/catallax-features/src/distance_build.rs), `positive`, `shell`, `thicken`; spline-shell test `unsupported_normal_envelopes_and_invalid_selection_keep_named_refusals` |
| Incomplete face/edge ownership proof or source topology coverage | Polygonal-loop, incident-face and source-name requirements can refuse; source code alone does not establish their observed frequency or that field ownership would resolve them | [planar shell](../crates/catallax-features/src/planar_shell.rs), `winding`, closed-loop/two-incident-face/named-source-face checks |
| Numeric proof/work refusal | `numbers.operation-budget` is a whole-request limit. `features.precision` preserves an undecided sign; planar-shell `distance.unsupported-envelope` may also carry an undecided classification | [sign adapter](../crates/catallax-features/src/profile.rs), `sign`; [planar shell](../crates/catallax-features/src/planar_shell.rs), `Outcome::NotExact` path. Keep diagnostic detail; ambiguous records stay unclassified |
| Downstream measurement coverage | A successful offset can subsequently refuse `measure.budget-exhausted` because finite curved trim integration is unsupported; this is not automatically an exhausted construction allowance | [public offset control](../crates/catallax-features/tests/trimmed_planar_offsets.rs), `a_spatial_reentrant_vertex_uses_a_spherical_distance_envelope` |
| Actual execution failure or incomplete record | `internal-error`, interrupted spans and missing responses must remain separate from conservative geometry refusals | [verification failure contract](VERIFY-INTERFACE.md#operation-log) and adjacent response-status rules |

## Follow-up priorities and coverage limit

1. **Obtain the actual current log bundle.** Required: public/authorized input
   identifiers and hashes, UTC run interval, kernel revision and binary identity,
   request allowances, per-adapter operation logs and final outcomes. Preserve
   missing/interrupted requests in the denominator. This is the blocking input,
   not a request to raise budgets or run a broad replacement campaign.
2. **Separate coverage from ownership and numeric limits.** The only observed
   selected events here are old unsupported-operation refusals. Current source
   now has shell construction routes, so these two records cannot support a
   present-day implementation priority. Preserve code, detail, operation stage
   and available undecided evidence; do not classify every `unsupported-envelope`
   as a surviving-face problem or every `budget-exhausted` as work exhaustion.
3. **Decide R11 step two only from a current classified share.** Report unique
   public cases as well as execution counts, with unsupported scope, invalid
   selection, proof budget, downstream measurement, actual failure and unknown
   records kept separate. There is presently no measured surviving-face share.

No kernel code, algorithm, allowance, tolerance, fixture, saved answer or kernel
output changes accompany this documentation. Historical observer/oracle results
embedded elsewhere in the checkpoint are not used as evidence in this census.
