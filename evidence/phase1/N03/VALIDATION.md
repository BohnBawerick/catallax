# N03 executed validation

Main/parent frozen: `ecad1812b339bb759d83ba8f057c927d69d393ec`.
Branch: `phase1/N03-approximate-step`. Date: 2026-11-15 UTC.
N02 merge `3e2a929e709c7916bfa8519e6cf5f7179cbc7cd5` and the supplied E01 merge
`ecad1812b339bb759d83ba8f057c927d69d393ec` were explicitly verified as ancestors.
No competing branch/PR was found before work; earlier work and E03 files remain.

## Portable reproduction

Run from the repository root using existing authorized Rust and OCP tools.
Set `N03_CAPTURE_DIR` to the absolute local location of `evidence/phase1/N03`
because Cargo tests run from their crate directory. The test writes the two
numerical files and evidence and checks their actual reread bytes.

```sh
CATALLAX_N03_CAPTURE="$N03_CAPTURE_DIR" cargo test --offline --release -p catallax-step --test numerical_step -- --nocapture
cargo test --offline --release -p catallax-step --doc numerical::NumericalStep
cargo test --offline --release -p catallax-step --test contract case01_block_fixed_entities_and_exact_counts -- --exact
python3 tools/step-regressions/check_numerical_step.py evidence/phase1/N03/two-loop.step evidence/phase1/N03/repaired.step --out evidence/phase1/N03/raw-reference.json
ci/run fmt
ci/run clippy -p catallax-step
git diff --check
```

Actual Rust commands selected existing authorized Cargo/Rustup roots and a reusable
Cargo target cache via local environment variables. No toolchain/profile,
permission, credential or global policy was changed. Default Python already had
OCP 7.8.1.1 available. No installation, alternative oracle or personal machine was
used. The unchanged dual-bin manifest warning remains diagnostic.

## Actual results

Six focused Rust tests pass. Full source and emitted payload numerical evidence
are captured in `two-loop-evidence.txt` and `repaired-evidence.txt` (actual Rust
debug records, not a stable wire schema). JSON files are explicit test summaries,
not exact/bounded result manifests. Source evidence retains the 0.01 mm local
tolerance and zero/one actual repair respectively; emitted repair histories are
empty. Scalar rounding limit is independently 1e-7 mm, measured maximum zero.
No complete body bound is declared. Both files have 11,325 bytes and export work
180,462 under a single 20M request allowance. Combined fixture/source-sew/export
work is 349,106; allowance-1 refuses. These are serial operation counts, not
worker comparisons or timing/performance qualification.

The two numerical STEP files were actually written and read from disk. The
artifact checked complete byte identity and planar seam syntax after reread.
An altered byte refuses `step.numerical-file-changed`. Other tests cover tolerated
source-gap refusal with unchanged input, decimal collapse refusal, unsupported
flat spline refusal (no plane/mesh substitution), invalid options and exhausted
allowance. The separate compile-fail conversion test passes.

The exact unit-block writer was captured as `legacy-block.step` before successful
numerical writer validation; existing exact-path source was never changed. Its
full 6,360 bytes remain identical in the final Rust test, at 2,400 work units and
under unchanged options/allowance. Required and legacy seam options produce the
same bytes. The independent existing exact STEP block contract control passes.
This is not the entire saved-answer or strict STEP replay campaign.

Official `ci/run fmt` passes formatting and the tooling suites: 88 checks tests,
15 edge-hunt tooling tests, 60 release-tooling tests. These use synthetic tooling
fixtures, not a geometry replay campaign. Scoped official Clippy passes all STEP
targets with warnings denied. No manually launched hosted CI is involved.

## Actual raw reference checks

`raw-reference.json` records OCP 7.8.1.1, Python 3.12.14, exact file hashes,
healing-disabled settings and complete two-case results. Both actual files import
raw as valid solids, with 16 vertices, 24 edges, 10 faces, two two-loop faces,
Euler-with-holes 0 and volume approximately 12 mm³. After the specified box cut,
validity and topology remain correct, volume is 10 mm³, the removed/retained wall
and hole probes have the expected material states, and the new boundary is found.
The files are unchanged by the reference check. Current file hashes were checked
again after final Rust captures and still match the reference records.

The first invocation of the older shared raw helper failed before geometry
validation because OCP 7.8 requires an explicit `StepData_Factors` argument.
The N03-only script uses that argument on the same direct raw builder with
identity millimetre/radian factors and the same 1e-7 precision/max tolerance.
It never calls TransferRoots, ShapeFix, SameParameter or a healing/sewing pass.
This was a binding-signature correction, not a tolerated geometry failure or a
switch of reference kernel. The repaired source control retains its N02 repair
history; raw refers to reference import without additional repair.

## Provenance, hygiene and limits

`artifact-digests.json` hashes current public source, script, documentation and
all retained artifacts, plus selected baseline/N02 fixture inputs. This manifest
is not self-hashed. Public prose and captures contain no machine paths. Full
numeric records are scoped evidence, not complete geometric bounds. Private
original tooling logs and reusable cache are retained; only disposable test
scratch is removed by tests. Earlier card evidence and worktrees are untouched.

No full P02, saved-answer or broad duplicate campaign was run. General curved
numerical export, arbitrary non-axis planes, merely tolerated source gaps, global
3D self-intersection/separation, native app round-trip/edit chains and all-worker
qualification remain outside this evidence. Exact/bounded producers, proof
boundaries, saved answers and global tolerances/allowances remain unchanged.
