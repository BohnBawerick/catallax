# Frozen admission export validity

The complete capture records **322 of the original 3,000 files (10.7%) with
STEP exports valid after OCCT builder internal pcurve healing**, beside **2,701 of 3,000 (90.0%) healed-valid**
results from the existing same-source OpenCASCADE reference.

The conditioning step is **400 admitted files, of which 322 are valid after builder healing
(80.5%), 61 invalid after that construction, 16 export-refused and one timed out**. All 400 imports
reproduced their expected body counts. The frozen inputs contain 808 bodies; one
valid body never qualifies an otherwise invalid multi-solid file. All ten
[retained P06 controls](retained-real-controls.json) are valid after builder internal healing.

[Results](results.csv), the [actual build/capture manifest](results.json),
[reference comparison](comparison.json), [authored controls](authored-controls.json)
and [separate byte/unit audit](independent-byte-audit.json) retain the evidence.
The 16 writer refusals are 15 `step.numerical-gap` and one aggregate export work
limit. The timeout is `00004777`, after import completed, during export; its
export work is unknown and it receives no validity credit. Raw invalidity is a
structural reader result, not a diagnosis of its underlying cause.

The observer producer is `4b422bcf65092239f0eeb4581d916595e7697560`, with production
kernel/exporter fixed at `f87c18ab1beb5e6fa74d74e9d0f06b8a24db0a38`. Its actual Cargo
artifact is SHA-256
`9f7648a9fe63060a274889d1f261fd5e80502c05c208af101134ced8299ecd77`.
The result table hash is
`7a8556150ba70f17d4bd15a7ca11fc7750822407b2469aaae22ab910f718ad69`.
No new admission rule, surface-body result or reference rerun enters this
comparison. The complete 3,000-file admission capture remains frozen separately.

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
60-second observation covers native import/export and the subsequent OCCT builder/reader
for the whole file. Each native request has one kernel worker. Builder reads run
in separate processes, avoiding shared OCCT reader settings between file workers.

`--resume` requires identical input/build/protocol identity. Existing terminal
observations, including timeouts, remain terminal. Process errors may be retried
into a new private attempt directory. Prior captures are preserved. An export
timeout cannot erase the flushed import checkpoint.
The source digest is checked when the producer exits, including either timeout
path.
Observer write, reread and byte-integrity faults are process errors and cannot
satisfy census completeness.
These collector safeguards were added after the frozen capture; its recorded
refusals contain no observer faults, and its data and producer identity remain
unchanged.

## What qualifies a file

The observer reuses
[the numerical STEP protocol](../../../tools/step-regressions/check_numerical_step.py):
explicit `StepToTopoDS_Builder`, millimetre/radian factors `(1,1,1)`, precision
and maximum tolerance `1e-7` mm, with no observer-added root transfer, ShapeFix,
sewing or SameParameter pass.
This does not disable processing inside the OCCT builder.
OCCT 8's [edge-loop construction](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateEdgeLoop.cxx#L787-L818)
projects vertex parameters through `ShapeFix_EdgeProjAux` and invokes pcurve checks.
Its [pcurve checker](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/DataExchange/TKXSBase/XSAlgo/XSAlgo_ShapeProcessor.cxx#L317-L468)
can use `ShapeFix_Edge::FixSameParameter` and rebuild or replace pcurves.
The reported validity therefore includes OCCT's internal healing; it does not
establish raw or unhealed validity of the emitted STEP geometry.
Each exported body must produce one explicit solid, and
`BRepCheck_Analyzer` must find its whole shape valid. Actual file hashes and
unchanged input bytes are checked independently of the observer's declaration.

A file qualifies only when the complete exported body set matches its frozen
admission and **every** body is valid after builder internal healing. One valid body cannot qualify a
multi-solid file.
The frozen capture's `raw_*` outcome/field names and `exported_step_healing:false`
and `exported_step_same_parameter:false` flags are preserved byte for byte.
Those legacy names/flags describe this direct-builder protocol and the observer's
lack of an additional repair pass; they do not assert that OCCT performed no
internal healing.
Outcomes remain distinct:

- `raw_valid`: all expected bodies were exported and pass the post-builder solid check.
- `raw_invalid`: an exported body fails the post-builder whole-shape validity check.
- `raw_reader_refused`: the builder/reader cannot construct the required solid.
- `export_refused`: a production writer returns a refusal, including work limits.
- `reimport_refused` or `admission_changed`: the frozen admission cannot be
  reproduced with its expected body set.
- `timeout`: the fixed whole-file observation cap expired; its stage is retained.
- `process_error`: observer I/O/integrity faults, missing results, protocol/schema failures or crashed processes
  prevent a complete publication.

The output table retains all 3,000 original IDs, their original first refusal
and frozen after-fix cause. The 2,600 inputs outside the full-admission population
retain their original partial/refused/timeout interpretation and receive
`outside_frozen_full_population` in the validity column. Partial admissions are
not promoted to full results.

`complete:true` requires all 400 terminal observations, matching source hashes,
the entire expected body set for every builder-valid file, and no process errors.
An incomplete result has no headline percentage. Public extracts omit source
geometry, raw diagnostics, private archive member names and workstation paths.

## Interpretation

The final comparison must publish the count valid after builder internal healing out of 3,000 beside the
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
