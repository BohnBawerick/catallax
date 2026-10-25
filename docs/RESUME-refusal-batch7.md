# RESUME: refusal repair batch 7

## Parked

Batch 7 remains **0/9 repaired** for the unfinished-coverage targets in
[PR34's remaining-nine table](https://github.com/BohnBawerick/catallax/pull/34).
[PR48](https://github.com/BohnBawerick/catallax/pull/48) stays draft. Stop
implementation and testing. Resume only on explicit instruction after the
Boolean-after-edit repair lands; there is no automatic restart.

Published PR48 head: `644c38834caf8f62a25b73c873ff728e3f0a20a1`.
Last inspected main: `3d59965621e3106997b45b5249a7be2f615e0b54`.
Their integration conflicts are in `catallax-numbers/src/algebraic.rs`,
`decider.rs` and `operation.rs` under `crates/`. This note does not resolve them
or qualify the published source against that main.

## Frozen baseline evidence

Baseline commit: `5b5477e433f89952e6d4806917fd5f273c835c3f`.
Its CLI binary SHA-256 is
`8c4791fbf4e5de12273fbf48de8565af92b78d180cb2d8d38e64e2a6534808c6`.
The verification driver SHA-256 is
`b0d0042dcf90f60d8a9304ccd3a75c03c36cad4d777013c7f43ffde27c578d53`.
Runtime source and Cargo dependencies are unchanged between the frozen base
and the last inspected main; later release tooling, documentation, site and
app-test changes are separate. Keep these evidence identities distinct.

Every row completed independently at **1, 2 and 4 workers**, with no cache hits:

| Approved input set | Request allowance | Answered | Refused |
| --- | ---: | ---: | ---: |
| Prior 477, historical completion envelope | 1,000,000,000,000 | 468 | 9 |
| Same 477, verification default | 100,000,000 | 423 | 54 |
| Supplemental 32, verification default | 100,000,000 | 28 | 4 |

Worker comparisons preserve strict geometric manifests, STEP bytes and artifact
presence. Complete responses, work and operation logs agree after excluding
only wall time. These counts distinguish answers from unchanged refusals; the
historical envelope is not default-budget acceptance. No allowance, precision
or tolerance was increased to claim a repair.

The nine targets all refuse `numbers.operation-budget` at the default. Seven
Boolean targets need 0.52–1.70 billion units to complete on the frozen baseline.
The two skew-cylinder measurements consume 209 million and 766 million units
in measurement alone; the latter still refuses `measure.centroid-not-separated`
at the historical envelope.

All eleven self-touching-shell controls were confirmed by baseline runs and a
rebuild of the published historical head. That historical binary SHA-256 is
`a1cca1630d7e0afba72704cb653fbb27296ea9e54cd2d04e8179843b01c0f349`.
This does not establish control preservation for an unqualified experiment.

## Five unsuccessful experiments

Each fresh default-budget probe answered **0/9**. Each retained all nine strict
refusal manifests and produced no STEP artifact. The archive commits below
are unpublished experiments, not proposed PR48 changes or qualification heads.

| Experiment | General change | Local archive commit |
| --- | --- | --- |
| 1 | Request-owned sign, enclosure and generator-reduction memos | `ed3ac7a0c0212a4707e570569730a1ddd98ccbb8` |
| 2 | Quotient-operation memos retaining local cold-proof fees | `2f388a781d4361d0dee7eb74cf347defa38c3dcb` |
| 3 | Exact generator and quadratic-tower identity shortcuts | `3045febd7552a677c0dddfa472dc878d9b0ef11e` |
| 4 | Deferred fractions with a final denominator-unit proof | `153fdcfea8bb66c8f3e21b7811ac1c92fafe6610` |
| 5 | Request-owned certificates for individual interval arithmetic steps | `dfd3e06094e32a666524d5e21e46472d58c6e076` |

Frozen probe binary SHA-256 identities, in the same order:

1. `38590cce8a09859e7985d8a7630e80ce9cd3aa1e717bfe55d8fb0a6fe0971152`
2. `aa4e5189b476d0340e596ab39162b119013dfe22e6a281eebe7696afb1b49343`
3. `f66a9f8bba418d4e09d9f4e39bf15edcacceaff36fe92338fe074a9dd72baf1f`
4. `d06d9b64f477d581ee238c7b585f6b060bad17701d095f4fedb32a3007470a8b`
5. `bc22d0db16d9129af1193249ba71fb8ebb96575a487b4f8b4d4fd3df1264aa40`

All experiments started from `68f0e66572d6986e0e5e0edab4b68bc0e0a5d8cb`.
Only experiment 4's CLI was built from its clean archive head. The other probe
binaries were built before archival, with recorded source patches; their
archive commits are not exact-head CLI validation. Experiment 4 passed five
focused tests. Experiment 5 passed all 60 numeric library tests, including five
new cache tests. No experiment completed a full workspace suite, full-477
candidate preservation matrix, candidate control matrix or hosted CI.

## Proposed repair map when explicitly resumed

- **Shapes:** retain exact carrier membership and endpoint-equality evidence
  through point construction, import and copying; inspect `store.rs` and its
  Boolean callers before changing ownership. Reuse definition-backed facts,
  never proximity or a small residual.
- **Boolean:** retain correlated meeting, germ/contact ordering and whole
  probe-path certificates through classification. Inspect `geometry.rs`,
  `probe_path.rs` and region/project callers. Coordinate with the Boolean edit
  repair before touching shared files. Scalar memos alone did not close this
  work. The [OCCT Boolean specification](https://occt3d.com/dev/doc/overview/html/specification__boolean_operations.html)
  retains incidence and interference records; that suggests an investigation,
  not permission to import its tolerance-based acceptance rules.
- **Measurement:** investigate reusable trace inputs and panel certificates
  across property refinements in `engine.rs` and `integrate.rs`. Numeric methods
  account for less than 5% of the skew targets' default-budget work. Preserve
  rounding, retry order, property budgets and published enclosures.

At restart, reconcile ownership and the intended main/source combination first.
Preserve all prior evidence. Require genuine target answers, unchanged earlier
answers and controls, strict manifests/STEP, 1/2/4-worker agreement, the full
suite and exact-head CI before claiming readiness. The old published head's
nine green checks do not qualify a new integration. Keep proof inputs and raw
evidence outside public source; no private inputs or evidence locations are
recorded here. No experiment source is included with this note.
