# Import correspondence: first kernel step

Producer: `a034ecf032e32511f5c58155eb945e18822ed59c`, based on main
`e9887bad203e351cd138b3eb8003e4202400cc05`. This is one step toward
ROADMAP 1.2. The [public contract](../../interface-changes/2026-11-15-import-correspondence.md)
describes the verifier, accepted derivative snapshots and unresolved cases.

## Real controls

The ten identities are the fixed P06/M5 answering controls, in their existing
order. They were recovered from the published checksum-pinned ABC chunk using
seed `20261111`; all 3,000 selected identities matched the historical census CSV.
All source files were individually hashed. This run evaluates ten correspondence
controls, not the 3,000-file import census.

Each replacement shifts every actual Part 21 entity identifier and reference by
1,000,000,000. Strings, comments, coordinates, parameters, orientations and other
geometry bytes remain unchanged. Both sources are independently admitted through
`apply` before comparison. Every verified target key was also checked against
the known renumbering. The standalone production example is
`cargo run --release -p catallax-import --example replacement_correspondence -- OLD.step NEW.step`.

The [machine-readable results](real-controls.json) retain old/new source hashes,
the actual binary hash, counts and charged work per stage. All 53 faces and 82
physical source edges match; no old or new entry is unresolved or unmatched.

| Public model ID | Verified faces | Verified edges | Comparison work |
| --- | ---: | ---: | ---: |
| 00000131 | 4 | 4 | 4,654 |
| 00000173 | 5 | 4 | 5,296 |
| 00000190 | 6 | 12 | 13,268 |
| 00000405 | 4 | 4 | 5,166 |
| 00000511 | 6 | 12 | 22,374 |
| 00000605 | 6 | 12 | 12,170 |
| 00000632 | 6 | 6 | 7,834 |
| 00000633 | 7 | 14 | 16,514 |
| 00000695 | 6 | 12 | 13,320 |
| 00000821 | 3 | 2 | 2,906 |

Each stage retains the existing 20M allowance. The table reports comparison
work; the JSON separately records both imports. It does not claim a combined
request was run under one allowance. These are fresh executions on a shared Mac
with two kernel threads, using a release binary rebuilt at the producer commit.
No latency qualification is claimed.

The source and unresolved per-model licence scope remain those of the
[original ABC census](../../milestones/R8-public-step-census.md#source-and-licence-scope).
Only public IDs, hashes, counts, outcomes and charged work are published. Actual
source files, renumbered files and full per-entity reports remain in private task
evidence. No external CAD oracle is needed for exact definition equality.

## Validation and remaining scope

Local fallback during the work-PC outage: all 259 import all-target tests and
20 feature import/replacement/selection tests passed with four build jobs and
four test workers. CI-selected all-target Clippy passed for the import crate
and its reverse consumers. Rust formatting passed. The existing local Python
format/tooling gate requires Linux affinity and GNU realpath, so its complete
result is provided by the restored work-PC run instead of being called a local
pass. The full Linux `ci/run fmt` gate and all 259 import all-target tests passed at
the producer commit on the work PC.

Controls cover changed geometry behind reused entity numbers, duplicates in
either source, reordered solids, actual planar derivatives and reopen, declined
solids, unavailable internal spline derivatives, public-record tampering and
exhausted comparison work. Existing replacement tests cover the editorial
transaction and rollback independently of this new report.

The report supplies verified matches and unresolved entries. It does not yet
wire automatic correspondence into `replace_step`, qualify the SolidDraft UI,
admit additional numerical spline families or establish a new import share.
Internally repaired spline results remain unresolved until their accepted finite
definitions can be compared. The full census and the app transaction retain
ROADMAP 1.2's remaining gates.
