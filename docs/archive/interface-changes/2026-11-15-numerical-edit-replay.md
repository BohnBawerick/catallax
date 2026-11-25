Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Numerical edit and replay history

Piece: numerical topology, feature evaluation and approximate STEP reports.

`NumericalBody::rebuild` checks a new definition under explicitly supplied current
controls. `replace_surface` and `repair_vertex` now keep prior checked states in
`history()`. Each entry owns its definition, tolerance and evidence. Current
`evidence()` contains current checks and cumulative repairs. Historical samples
never become current checks or a complete bound. `source_definition()` keeps the
original input, including endpoints removed by explicit sewing.

`catallax_features::evaluate_numerical` evaluates a `NumericalRecipe` separately
from exact recipes and checkpoints. The recipe retains an immutable source,
explicit sewing pairs, source tolerance and ordered edits with their own controls.
The first edits support coordinate-axis dimension changes on rational planar
bodies with unit coordinate frames, and half-space cuts across two-layer straight
prisms. Cut admission checks retained line directions and exact axial endpoint
coordinates as well as nominal vertex alignment. Tilted or displaced traces refuse.
Cuts move the selected cap, retain its source origin and add an edit cap
origin. The finite outer and hole loops remain explicit. Unsupported geometry,
nonpositive dimensions and out-of-range cuts refuse. This is a bounded operation
family, not general numerical Boolean coverage.

`NumericalRecipe::to_json` and `parse` save generating definitions as
`catallax-numerical-recipe/1`. `replay_numerical` parses and rechecks the source,
sewing and every edit under one active request allowance, or the unchanged 20M
default. Input bytes are charged before JSON decoding starts. Saved checks are not trusted. `NumericalBuilt::edit` borrows the previous
result and returns a complete new result only on success. Failed checks, precision
limits and exhausted work return no partial result.

`NumericalStepReport` adds `original_definition` and `source_history`. Its existing
source evidence still describes the current numerical body; emitted evidence
still describes the separately checked decimal file. Existing exact APIs and
exporters need no migration. No tolerance, precision cap or work limit increases.

The opt-in driver is:

```sh
cargo run --release -p catallax-features --example numerical_replay -- saved.json result.step
```

It shares 20M work across replay and STEP construction. It publishes approximate
STEP only after those checks pass. It writes a unique adjacent staging file,
rereads and verifies its bytes, then atomically renames it to the output path.
Input/output identity is checked before publication and refusal cleanup, including
symlink and Unix hard-link aliases. Refusal removes stale output only after that
guard passes. Staging files are removed on success and failure. The JSON summary reports current tolerance, prior-stage
controls, repair counts, origins and deterministic work. This example does not add
a numerical route to the existing request protocol or the app checkpoint system.

Regression commands:

```sh
cargo test --release --locked -p catallax-topology --test numerical_trim
cargo test --release --locked -p catallax-features --test numerical_edits
python3 tools/step-regressions/check_numerical_edit.py edited.step
```

The last command needs the existing OCCT Python dependency. It uses the raw
healing-off builder, checks volume and material on both sides of the new cap,
checks the retained through-hole, then cuts the imported solid once more.
