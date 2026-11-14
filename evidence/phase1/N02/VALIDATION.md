# N02 executed validation and retained evidence

Baseline: `f26c898c94582dd640c754b210cca796071416d4`.
Branch: `dot/n02-numerical-trim-sewing`.
Date: 2026-11-14 UTC.

Portable command forms (run from the repository root):

```sh
cargo test --offline --release -p catallax-topology
cargo test --offline --release -p catallax-app --test numerical
ci/run fmt
ci/run clippy -p catallax-topology
```

Actual commands selected existing authorized Cargo/Rustup roots and a reusable
Cargo target cache through environment variables; they did not change the
repository's build profiles or toolchain. `cargo` was invoked from that existing
tool root. No toolchain installation, permission change or credential action
occurred.

Final scoped capture, with `CATALLAX_N02_EVIDENCE` set to the absolute local
location of this repository-relative directory:

```sh
CATALLAX_N02_EVIDENCE="$N02_EVIDENCE_DIR" cargo test --offline --release -p catallax-topology --test numerical_trim -- --nocapture
```

Observed final focused run: exit 0, three tests passed. Relevant
transcribed output:

```text
N02 control=crossing refusal=numerical.invalid-trim
N02 control=touching refusal=numerical.invalid-trim
N02 control=outside refusal=numerical.invalid-trim
N02 control=reversed refusal=numerical.invalid-trim
N02 control=self-crossing refusal=numerical.invalid-trim
N02 gap=0 result=numerical work=168612 exhausted=false repairs=0
N02 gap=0.005 result=numerical work=168612 exhausted=false repairs=1
```

Both successful work records come from `with_operation_budget(20_000_000, ...)`
around sewing and the entire admission, on the executing thread. They are not
worker-count comparisons. The separate allowance-1 control exhausts and refuses.
No global allowance, work charge or tolerance was increased. Test scheduling may
run independent test functions concurrently; the operation scopes themselves are
serial. There is no production parallelism or multiworker qualification claim.

The complete topology suite passes 59 tests plus one compile-fail documentation
test. The existing N01 application numerical suite passes all 10 tests unchanged.
Official `ci/run fmt` passes workspace formatting and 88 tooling-check tests,
15 edge-hunt tooling tests and 60 release-tooling tests. These tooling tests use
their own temporary synthetic fixtures; they are not hosted gates or a kernel
replay campaign. Official scoped `ci/run clippy -p catallax-topology` passes
`--all-targets -- -D warnings`. The existing workspace dual-bin manifest warning
is diagnostic and unchanged.

The retained `near-gap.txt` is actual final `NumericalDefinition` and
`NumericalEvidence` Rust debug output. It includes the two loops, every required
sample/enclosure, the stated local tolerance, completed checks, removed-to-keeper
mapping and one actual vertex-movement record. The other text records are actual
refusals. They contain no partial body. These debug records are not a new wire
schema and grant no exact, bounded or export status. Source/record hashes are
listed in `artifact-digests.json`.

One initial optional capture attempt failed because its repository-relative
output directory was resolved from Cargo's crate working directory. Supplying the
existing absolute local destination corrected the evidence sink; no geometry,
acceptance rule or production code changed for that correction. The recorded
final capture and tests pass. Public records use repository-relative references
and contain no machine paths. Original tooling logs are retained privately;
this file is an explicit summary, not a purported raw stdout log.

Static checks: relative document links, published text machine-path audit,
changed-artifact digests and `git diff --check`. No new build/scratch directory
was created beyond the reusable Cargo cache and the required private tooling-log
capture. User checkout/work, previous branches and required evidence remain.

Limits: no rerun of P02's M10/M8 requests, no 477-answer replay campaign, no manual
CI, no C++ oracle or private input, no new numerical STEP file, no independent
import or full application edit chain. N01's existing sampled-check limits and
generic-constructor trust caveat remain. Exact material/volume queries and global
3D self-intersection/separation are not established by this planar family.
