Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Verified STEP file transaction

`catallax_features::import_replace` adds:

```rust,ignore
replace_step_file(
    recipe: &Recipe,
    edit: &StepFileReplacement,
    files: &ImportFiles<'_>,
    dec: &mut Decider,
) -> Result<FileReplacementOutcome, Refusal>
```

`StepFileReplacement` supplies the expected complete old hash, a new `FileRef`
whose immutable bytes are already in `ImportFiles`, and new hash-bound import
decisions by affected step. All STEP import steps using that old hash participate
in the transaction. Existing approvals are not copied to different bytes.

`Ready(VerifiedFileReplacement)` returns a fully rebuilt revision and analysis.
`NeedsResolution(FileReplacementAnalysis)` returns correspondence and unresolved
JSON-pointer occurrences without a candidate model. Unverified faces, source
solid membership and wildcard/non-imported origins are explicit unresolved kinds.
Other failures are typed refusals. No failure changes the borrowed old revision.

The analysis includes a hash of the canonical old recipe. Before swapping its
document, a caller verifies that revision is still current. Save the returned
complete recipe and retain immutable old/new bytes for reopen and undo. The app
owns document commit and explicit unresolved choices; reports are not authority
accepted back by the kernel.

The existing `import_edit::replace_step` contract is unchanged. The shared
hash-bound decision decoder is extracted from recipe evaluation without changing
its accepted decision format. Each rebuild meters source bytes before hashing,
once per distinct digest, including unchanged imports and retained prefix sources.
Source and rebuild refusals keep their underlying codes and undecided evidence,
with source-stage or failing-step context in their detail. Exhausted requests
discard successful outcomes while retaining evidence on an existing refusal.
The new transaction targets existing STEP recipe
imports; it does not add a numerical-body cast into the exact evaluator.

See [qualification and reproduction](../phase1/import-file-replacement/README.md)
for ten real sources, complete source hashes and multi-solid/downstream controls.
