Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Verified STEP file replacement

`catallax_features::import_replace::replace_step_file` replaces every STEP import
step using one old source hash as a single candidate revision. It recomputes
[correspondence](../import-correspondence/README.md) from the actual old and new
accepted source definitions. Entity addresses, source-solid order, labels and
proximity do not establish identity.

A selected source solid follows complete, unique verified face membership.
Every structured origin in the canonical recipe is visited, including nested
copy/offset origins, face references, edge and vertex names, touching filters
and fillet setback references. Parsing the rewritten recipe restores canonical
name ordering; rebuilding the complete recipe checks downstream selections.
Body step identities and the order of directed edge selectors are preserved.

`Ready` contains the complete new build and recipe. `NeedsResolution` contains
correspondence and every unresolved structured occurrence, without a candidate
build. Missing sources, stale hashes, incompatible geometry and failed rebuilds
return typed refusals. The borrowed old recipe and source bytes remain intact.
An application compares `analysis.old_recipe_sha256` with its current revision
before committing the result and retains the old revision and bytes for undo.
Serialized reports are informational and are never accepted as transaction
authority. The complete transaction shares the unchanged 20M work allowance.
Rebuild source verification meters bytes before hashing and verifies each
distinct source once per call, including unchanged imports and cached prefixes.
Source and rebuild failures retain their underlying refusal code and typed
undecided evidence with stage context.

This is a kernel transaction for existing recipe imports. SolidDraft's document
swap, unresolved-reference interaction, history UI and numerical-body group
integration remain in the app slice. The earlier explicit editorial
`replace_step` API remains available with its existing source-selection rules.

## Actual controls

Producer: `27d0d8157b460f6f388ffdef74f9c92e91920cef` on the Linux timing system,
revalidated after metering all rebuild source verification and retaining typed
refusal evidence.
[real-controls.json](../../../phase1/import-file-replacement/real-controls.json) pins the executable and both source
hashes for all ten fixed P06/M5 capability parts. Only actual Part 21 entity
addresses were shifted by 1,000,000,000; geometry, strings and comments were
preserved. Both hashes match the preceding correspondence qualification.

All ten complete transactions passed. They checked 53 actual imported face
origins, 82 edge names and 40 vertex names against the rebuilt bodies. No model
edge or vertex was left unnamed in these controls, and there were no unresolved
references. Each saved replacement recipe reopens and rebuilds; the retained old
recipe also replays with its old origins. Each complete qualification, including
the original build and both replays, stays under the unchanged 20M allowance.
Measured replacement work ranges from 68,300 to 521,231 units; the largest
complete qualification uses 1,059,706 units.

Eleven authored controls cover the missing legacy selected-solid transaction,
two selected solids in one changed file, reordered source solids, saved/reopened
revisions, nested copy/edge/vertex/touching/setback references, fixed face and
adjacency filters, changed geometry behind reused numbers, duplicate geometry,
stale decisions, corrupt bytes and exhaustion, including a later unchanged
source with insufficient remaining work to hash its bytes. All 35 scoped release
import tests pass, including the legacy import, selector, IGES and seamed-file
controls. Three additional release unit controls retain underlying refusal codes,
typed undecided evidence and source/rebuild stage context, including exhaustion.
The Linux format gate and CI-selected eighteen-crate all-target Clippy pass.

Build the `catallax-features` example `import_file_replacement` in release mode.
With the unchanged private `MODEL_ID.step` sources, capture outside the checkout:

```sh
python3 -B tools/research/check_import_replacement.py --binary TARGET/release/examples/import_file_replacement --sources PRIVATE_SOURCES --out NEW_PRIVATE_CAPTURE
```

The runner verifies original and changed-file hashes before invoking the
transaction. Raw source copies stay private under the
[existing licence scope](../../milestones/R8-public-step-census.md#source-and-licence-scope).
The committed table contains only IDs, hashes, counts, work and outcomes.

The relation follows OCCT's
[`BRepTools_ReShape`](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V8_0_0/src/ModelingData/TKBRep/BRepTools/BRepTools_ReShape.cxx)
approach: an explicit verified replacement relation supplies substitutions.
Canonical complete source definitions establish the relation before the rewrite;
no nearest-face rule is used. The [P04 SMLib map](../P04/PORT-CONSUMERS.md)
requires no new spline projection or tracing port for this transaction.
