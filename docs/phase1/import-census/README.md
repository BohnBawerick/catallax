# Fixed-input numerical import census

The sample is the original 3,000 ABC model IDs, seed `20261111`, from
`abc_0000_step_v00.7z`. The archive SHA-256 is
`52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe`.
The [individual source hashes](source-hashes.json) have SHA-256
`3323e104c377023846073a05d7e164ca1082dc3dcf32cc046eee897401454cc7`.
They were recovered from the pinned archive and independently checked after
transfer to the work PC. Source geometry and archive member names stay private
under the [existing licence scope](../../milestones/R8-public-step-census.md#source-and-licence-scope).

The [published R8 baseline](../../milestones/R8-public-step-census.md) has
186 full source-coverage admissions (6.20%), three admissions with unverified
source coverage, 2,782 refusals and 29 external observation timeouts.
New observations use the opt-in numerical import API while preserving the
20,000,000 work allowance and the external 60-second observation cap.

`full` means the import returned bodies and the source traversal established
complete coverage. `partial` means bodies were returned with source coverage
unverified. `refused` means a typed import refusal. `timeout` means the external
observation cap expired, even if a checkpoint was written before termination.
Process errors and missing results prevent a complete report.

The numerical category remains approximate, including zero sampled residuals.
The importer requests the fixed 0.01 mm tolerance and records its actual local
source limits. It does not establish a complete-body approximation bound.
This census measures admission. The separate
[oblique STEP checks](../import-oblique-step/README.md) qualify actual pose and
export workflows; admission alone does not assert those operations for every file.

## Reproduction

Extract the original selected archive members into a private `SOURCES` directory.
Retain the private `[model_id, archive_member]` selection that reproduces the
original seeded IDs. Build `catallax-import`'s `r8_census` example in release mode,
then capture outside the checkout:

```sh
python3 -B tools/research/import_census.py --binary TARGET/release/examples/r8_census --sources SOURCES --selection selection-private.json --source-hashes docs/phase1/import-census/source-hashes.json --out NEW_PRIVATE_CAPTURE --workers 4
```

The collector verifies all 3,000 IDs against the original per-model table and
the pinned source-hash manifest. Each source hash is checked before and after
its observation. Four independent observers each use one kernel worker. Raw
source and import checkpoints remain private; the CSV and summary contain IDs,
hashes, counts, tolerances, work and outcomes.

A resumable capture binds the binary, observer, collector, source revision,
toolchain, source hashes, original table, selection, workers and caps. Use
`--resume` with that unchanged identity. A subset cannot produce a complete
report. The final table retains each original first refusal beside its new cause,
with transitions and recovered/lost full-admission IDs calculated from actual
observations.
