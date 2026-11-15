# B03 bounded inventory validation

Baseline: `48a6d53d57300d7b248e45a3aa86ebca6ad7f06d` (merged B02 #607).
Branch: `phase1/B03-general-drop`. Date: 2026-11-15 UTC.

Executed in Rust on only the two unchanged shipped `drefuse` rows:

```sh
cargo test --offline -p catallax-trace --test b03_drop_inventory -- --nocapture
ci/run clippy -p catallax-trace
rustfmt --edition 2024 --check crates/catallax-trace/tests/b03_drop_inventory.rs
git diff --check
```

Actual commands select the existing authorized Cargo/Rustup roots and reusable
build cache through local environment variables. The final test exits 0: one
test, two fixture cases, both integer refusal expectations pass. No positive
binary64 outputs are compared. Official scoped Clippy passes all targets with
warnings denied; formatting and diff checks pass. The unchanged dual-bin
workspace warning remains diagnostic.

Transcribed complete observations:

```text
B03 OFF direct general drop keep_all=false: curves=0 max_drop=0 max_deviation=0
B03 diagnostic reached DropTracer::trace_curve2 curve_count=0
B03 on_surface=true complete refusal=Some(NotPorted(Intersect))
B03 on_surface=false complete refusal=Some(NotPorted(Intersect))
```

The diagnostic line came from temporary Rust-only instrumentation of that
private helper. The source file was restored unconditionally and checked equal
to baseline; the final stock-source run reproduces both public refusals and the
empty direct-general result, without the diagnostic line. Original diagnostic
stdout is retained privately. The text above is an explicit transcription, not a
purported raw stdout file. No instrumentation is committed. Empty general-drop
results and zero-valued diagnostic maxima do not prove geometric emptiness or a
zero global bound.

Static inspection reads the 14 public P02 question records and both historical
worker runs. The structured inventory matches recipe names from portable
questions, retains complete refusing-operation records, and hashes references to
the complete source responses/operations. It does not run any P02 recipe and
does not add current branch reachability to those historical records.

No current Phase 1 consumer was dynamically mapped to OFF. No positive donor
output is public in the refusal row; no donor oracle was executed. Hence no
production fallback port or finite candidate curve acceptance is claimed. B02's
native-parameter guard and general-drop refusals are unchanged. No STEP/export,
application proof, worker comparison, allowance change or hosted CI run occurred.
Private diagnostic logs and reusable cache are retained as evidence; temporary
source edits were removed. Prior work remains untouched. Public paths/links and
artifact hashes are statically audited before publication.
