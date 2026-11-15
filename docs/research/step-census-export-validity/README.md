# Frozen admission export validity

This protocol qualifies the exported STEP files of the **400 admission-only**
results in the [complete periodic census](../step-census-periodic/results.json).
Those fixed inputs contain 808 bodies. These are population counts, not an
export-validity result. A headline requires a complete capture.

The original seed is 20261111 and the denominator remains 3,000. Input IDs and
individual byte hashes come from the original
[source manifest](../step-census/source-hashes.json). The frozen admission table
has SHA-256
`b08224a6f919405de46d35896d8d682febbd10c36495d34759f4958c80be3bd7`.
No new admission rule or surface-body result enters this comparison.

## Build and capture

The observer calls the unchanged production import entry at its default
numerical options, checkpoints its result, then exports every returned body
through the appropriate production STEP writer. Certified imports use the
existing model writer; numerical imports keep the approximate writer and its
evidence. No numerical body is cast into an exact model. Files are written and
read back before their actual bytes and hashes are recorded.

The builder requires a clean checkout and binds the actual executable reported
by Cargo to its source revision/tree, compiler, observer, collector and lockfile.
Every production source and dependency must match the supplied final reviewed
exporter revision; the additional example and research tooling may differ.
This keeps admission and exporter algorithms fixed during the measurement.
Build records and full Cargo logs stay outside the checkout. The public nested
build record is a whitelisted subset; `build_record_sha256` identifies the
original private record, not that subset. The measured producer remains distinct
from any later publication/rebase head:

```sh
python tools/research/export_validity_census.py build \
  --kernel-revision FINAL_REVIEWED_EXPORTER_HEAD --record PRIVATE/build.json
```

Run collection in a Python environment containing **OCP 8.0.1.1**, with the
original private selection and durable source corpus:

```sh
python tools/research/export_validity_census.py collect \
  --build-record PRIVATE/build.json \
  --sources DURABLE_SOURCE_CORPUS --selection ORIGINAL_PRIVATE_SELECTION.json \
  --source-hashes docs/research/step-census/source-hashes.json \
  --out PRIVATE/capture --workers 2 --timeout 60
```

Import has the unchanged 20M kernel allowance. The aggregate export of all
bodies has a separate 20M allowance. Both stage counters are retained; observer
I/O, metadata hashing and OCCT work are outside the kernel counter. The external
60-second observation covers native import/export and the subsequent raw reader
for the whole file. Each native request has one kernel worker. Raw readers run
in separate processes, avoiding shared OCCT reader settings between file workers.

`--resume` requires identical input/build/protocol identity. Existing terminal
observations, including timeouts, remain terminal. Process errors may be retried
into a new private attempt directory. Prior captures are preserved. An export
timeout cannot erase the flushed import checkpoint.

## What qualifies a file

The raw reader reuses
[the numerical STEP protocol](../../../tools/step-regressions/check_numerical_step.py):
explicit `StepToTopoDS_Builder`, millimetre/radian factors `(1,1,1)`, precision
and maximum tolerance `1e-7` mm, and no root transfer, ShapeFix, sewing or
SameParameter call. Each exported body must produce one explicit solid, and
`BRepCheck_Analyzer` must find its whole shape valid. Actual file hashes and
unchanged input bytes are checked independently of the observer's declaration.

A file qualifies only when the complete exported body set matches its frozen
admission and **every** body is raw-valid. One valid body cannot qualify a
multi-solid file. Outcomes remain distinct:

- `raw_valid`: all expected bodies were exported and pass the raw solid check.
- `raw_invalid`: an exported body fails the raw whole-shape validity check.
- `raw_reader_refused`: the raw reader cannot construct the required solid.
- `export_refused`: a production writer returns a refusal, including work limits.
- `reimport_refused` or `admission_changed`: the frozen admission cannot be
  reproduced with its expected body set.
- `timeout`: the fixed whole-file observation cap expired; its stage is retained.
- `process_error`: missing results, protocol/schema failures or crashed processes
  prevent a complete publication.

The output table retains all 3,000 original IDs, their original first refusal
and frozen after-fix cause. The 2,600 inputs outside the full-admission population
retain their original partial/refused/timeout interpretation and receive
`outside_frozen_full_population` in the validity column. Partial admissions are
not promoted to full results.

`complete:true` requires all 400 terminal observations, matching source hashes,
the entire expected body set for every raw-valid file, and no process errors.
An incomplete result has no headline percentage. Public extracts omit source
geometry, raw diagnostics, private archive member names and workstation paths.

## Interpretation

The final comparison must publish the raw-valid count out of 3,000 beside the
400 admission-only count and the existing **2,701 OCCT healed-valid** results.
The reference capture uses normal transfer/default healing and its own default
tolerances; it is a separate method, not a common-tolerance benchmark. Reuse its
complete same-source capture after auditing IDs and hashes rather than rerunning
or changing the reference protocol.

Raw B-rep validity does not establish source fidelity, a complete deviation bound,
material agreement, an arbitrary edit or a Boolean for every file. Numerical
results remain approximate. Those claims require their own evidence. The named
multi-seam import/pose/material qualification is retained
[separately](../../phase1/import-periodic-numerical/README.md).
