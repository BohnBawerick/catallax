## 2026-10-24: bounded export-only reconnection of collapsed planar wires

Pieces: 7, CLI verification

Change: `StepReport` adds `regularizations: Vec<StepRegularization>`. Each record
identifies the face and edge pair, the parallel neighbouring faces merged in the
file, and an exact whole-ribbon displacement bound. A contained straight wire
segment that overlaps another wire's segment in the actual decimal export is
reconnected only after proving a parallel exact gap within both the unchanged
export tolerance and 1 nm. The exporter uses an isolated working topology;
the source model, names, topology and manifest formats do not change.

Why: a positive exact separation can be smaller than the export's numeric
resolution. Separate inner and outer wires then coincide in STEP and fail native
topology checks. More decimal digits do not make that topology usable by a
binary64 reader.

Migration: `StepReport` struct literals must initialize `regularizations`, or use
`..StepReport::default()`. Counts in the report describe exported topology rather
than necessarily matching model counts. Both CLI export paths retain the bound
records in an optional `step-regularizations.json` sidecar; the verification
scoreboard reads its exported face count for OpenCASCADE checks. Existing
manifests and seam/ring counters retain their formats. The native checker does
not repair shapes. Valid, separated STEP wires retain their previous bytes.
Model-side gap merging remains separate tolerance-core work under PLAN item 11.

