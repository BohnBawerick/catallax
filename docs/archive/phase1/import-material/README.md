Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Material edit qualification after real oblique import

The fixed census part `00000022` numerically imports with the same source and
repaired-definition hashes as the [oblique writer qualification](../import-oblique-step/README.md).
After its history-preserving 3/5/7 mm pose, the existing N04 resize operation was
asked to double its X extent. That affine edit would double volume if supported.

The actual result is the typed refusal `features.numerical-unsupported`.
N04's current edit family requires rational lines and unit coordinate-axis plane
frames; this part has oblique plane frames. No material change was applied.
The protected geometry, original source definition, complete evidence and prior
history remain unchanged. Both the pre-attempt and preserved approximate STEP
files are byte-identical to the preceding writer's posed artifact. Each imports
raw through OCCT as a valid solid with unchanged topology, volume and interior
material. Its prior history stage is retained in both body and STEP report.

The paired supported control uses the authored N04 square tube. Doubling its
height from 1 to 2 mm changes material volume from 12 to 24 mm³, a known delta of
12 mm³. N04 retains the prior history, and saving and replaying the complete
edited numerical recipe reproduces its source, geometry, history and evidence.
Both actual STEP files import raw as valid solids. The reference measures the
volume delta, checks newly occupied material and confirms the hole remains open.

Producer: `65db0c50bc3b2988567215b153cd8769570dde2a` on the Linux timing system.
[controls.json](../../../phase1/import-material/controls.json) pins source, repaired-definition, executable and
STEP hashes, actual outcomes and work. The complete real attempt, preserved
exports, supported edit, exports and replay share the unchanged 20M allowance:
3,843,447 work units. Tolerance stays 0.01 mm; all results and files remain
numerical approximate with no complete-body approximation bound.

The reference uses OCCT Python 8.0.1.1 and the existing raw
`StepToTopoDS_Builder`, with explicit millimetre/radian identity factors,
precision and maximum tolerance of 1e-7 mm. It does not heal, sew, run
`SameParameter`, rewrite inputs or use `TransferRoots`. Authored volume checks
use 1e-9 mm³ comparisons. Real coordinates, requested absolute dimension and
geometric reference values stay private under the
[existing licence scope](../../milestones/R8-public-step-census.md#source-and-licence-scope).

## Reproduction

Build the `catallax-features` example `import_material_edit` in release mode,
then run the unchanged hash-fixed source into a new private directory:

```sh
TARGET/release/examples/import_material_edit SOURCE_00000022.step NEW_PRIVATE_RESULTS
python3 -B tools/step-regressions/check_import_material.py --results NEW_PRIVATE_RESULTS
```

The example checks source and repaired-definition hashes before the edit and
requires the pre-attempt STEP hash to match the preceding real writer check.
The reference independently checks actual STEP bytes, refusal preservation,
history metadata, material probes and the supported control's volume delta.
The full Linux format gate, both release example controls and the CI-selected
eighteen-crate all-target Clippy gate pass. This step qualifies the existing edit
boundary; it adds no modelling algorithm or writer family.
