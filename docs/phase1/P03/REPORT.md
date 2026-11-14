# P03 — remaining everyday feature baseline

Executed 2026-11-14 on public Catallax main **f4db345f2e7c51835d1ebeba9f6446d117144117**. Six requested cases give typed construction refusals. Five nearby controls answer and produce STEP files. No request exhausts its allowance. This is baseline CLI evidence, not feature qualification or an algorithm change.

## Provenance and method

Repository: `https://github.com/BohnBawerick/catallax.git`. `git ls-remote origin refs/heads/main` returned the baseline above before work and again after execution. A fetch and isolated worktree froze that exact commit on `phase1/P03-everyday-evidence` in an isolated checkout; the original checkout and older branches were preserved.

Read README, BUILDING, RECIPE sections 4.5.1/4.6/4.7, VERIFY-INTERFACE, relevant feature tests, fillet coverage and interface fragments. `.agents/skills/verify-catallax` resolves to `.claude/skills/verify-catallax`; its SKILL.md, feature map and answers instructions guided real CLI execution. No AGENTS.md was found in the repository. “Fillet/Wall/Bevel” separate decision documents were not present among visible repository sources; the public recipe definitions and cited implementations/tests below supply the available intent. No unseen decision was inferred.

The prescribed build command `nice -n 10 cargo build --release -p catallax-cli --bins -j 2` exited 127 (`cargo` not found). Used the already saved `target/release/catallax-verify`: doctor passed, `about.request = catallax-request/1`, embedded kernel commit equals frozen main, SHA-256 **1dfe54e616aee1aa353b3f6a07933ab4451f44a8791cf53ec4c748ffd6f48128**. A fresh rebuild and independent binary/source attestation remain unavailable; embedded identity alone is not that attestation. No installation or environment change was attempted.

Ran one fresh persistent worker pinned to allowed core 0, one kernel thread, ordinary recipe source, **20,000,000 work per complete request**. This is the existing application allowance, explicitly documented in VERIFY-INTERFACE; the driver verification default is 100M. No tolerance/export option, allowance policy, expected answer or feature algorithm was changed. The run returned exit 0; worker 0 exited 0; all 11 cases completed without cache reuse. Whole-request work includes export for successful controls. Timing is diagnostic on a shared machine, not a latency qualification.

## Executed outcomes

Exact rational inputs are in [inputs](../../../evidence/phase1/P03/inputs); [summary.json](../../../evidence/phase1/P03/summary.json) maps filenames to raw case keys, complete refusal details, work and topology. Each raw case retains request, response, manifest, operation log, portable question and artifact hashes; successful cases also retain actual STEP bytes.

| Requested input | Refusal at step | Code | Whole-request work | STEP |
|---|---|---|---:|---|
| `mixed-unselected.json` | round | `fillet.unsupported-corner` | 43,067 | absent |
| `partial-oblique.json` | round | `fillet.unsupported-end` | 485,948 | absent |
| `bored-side-opening.json` | shell | `distance.unsupported-envelope` | 13,180 | absent |
| `hollow-oblique-draft.json` | drafted | `draft.unsupported-oblique-annulus` | 454 | absent |
| `bevel-three-edge-corner.json` | bevel | `chamfer.unsupported-topology` | 224 | absent |
| `skew-hole.json` | drilled | `intersect.trace-native-cover` | 60,085 | absent |

| Existing-success control | Outcome | Faces / edges / vertices | Whole-request work | STEP |
|---|---|---|---:|---|
| `mixed-all-control.json` | answered | 38 / 74 / 38 | 724,536 | present |
| `bored-cap-control.json` | answered | 14 / 28 / 16 | 96,543 | present |
| `oblique-solid-control.json` | answered | 3 / 2 / 0 | 3,801 | present |
| `bevel-single-control.json` | answered | 7 / 15 / 10 | 2,695 | present |
| `radial-hole-control.json` | answered | 7 / 8 / 0 | 15,528 | present |

Refusals were returned by the requested construction operation after successful preceding steps, not inferred from missing files or a killed process. All six contain explicit construction coverage reasons, not impossible-geometry verdicts. The five controls' STEP existence and `ISO-10303-21;` headers were checked against actual bytes. Consult raw manifests for all topology fields; a file header is not an independent STEP validity proof.

## Exact case derivations and source interpretation

- **Unselected mixed vertex:** `parts/fillets/cases.jsonl`, `l-mixed-all`, reproduced first unchanged. The refusal input changes only its fillet selection to base/end–wall/y_max and base/end–base/x_min, exactly as `fillet_mixed_vertices.rs::a_mixed_vertex_with_an_unselected_edge_keeps_its_refusal`. Dimensions are base 20×20×5 plus wall 20×5×10, radius 1. Observed detail: “a convex and a concave strip meeting at a vertex need a general corner patch; maximum unknown”. RECIPE 4.5.1 and the 2026-11-08 mixed-vertices fragment distinguish full mixed selection success from partial selection coverage.
- **Partial oblique end:** unchanged public `partial-oblique` recipe, radius-5, height-16 cylinder, inclined cut at z=8 (quaternion 24/25,0,7/25,0), intersected with the x≥0 clip block, fillet radius 1/2 on cylinder-wall/cut-start. Observed detail: “vertex trimming needs straight rails; maximum unknown”. The coverage note identifies clipped ellipse ends on the vertical wall as missing; a closed oblique ring would not reproduce them. Public row metadata/oracle expectations were not executed or used as proof. Work differs from older source-reported 512030 because this report records the current executed binary/route, not the old result.
- **Bored side opening:** unchanged `fixtures/bored-shell.json` cap-open control: stock 40×30×12, diameter-6 through hole at (20,15,12), shell thickness 2. The requested case selects only stock/x_min as open instead of stock/end. This explicit one-face side selection specializes the side-opening family described in `bored_shell.rs` and the 2026-11-02 separated-through-bore-shell fragment; it is not copied from the broader `{step:stock}` regression selection. Observed detail: “bored shells open only their perpendicular end caps; maximum unknown”.
- **Hollow oblique draft:** `fixtures/draft-oblique.json` unchanged control, then append radius-2 bore circle at the origin to its radius-8 sketch. Exactly mirrors `draft_oblique.rs::hollow_oblique_refusals_name_the_unproved_gap_with_source_numbers` with only outer wall selected. Extrusion height 2; pull (0,0,1); neutral point (0,0,1), normal (1,2,8); angle 45°. Observed missing whole-height/whole-angle radial separation and shared trimmed cap rings; radii squared 64/4, cap span 2, tan angle 1, transverse ratio squared 5/64, certified radial gap unknown. No invented gap bound.
- **Three-edge bevel corner:** public `chamfers.rs::block_recipe` stock 10×12×8, distance 2, ordinary equal-distance roll-over. Single x_min/y_min edge is the control. Requested case adds x_min/start and y_min/start, the three edges at (0,0,0). This explicit extension uses RECIPE 4.6's existing corner definition, not a saved proof answer. Observed detail: “three selected strips at a vertex need a corner patch”. Source `chamfer_build.rs` contains the same coverage exit. A defined corner does not imply an implemented corner.
- **Skew hole:** exact unmoved `hole_neighbours.rs::collar` case from `a_hole_meeting_a_bevel_keeps_its_typed_refusal`: revolve rectangle with corners (6,0),(12,0),(12,12),(6,12) about y through 360°, chamfer outer bottom/top rims distance 1, drill diameter 2 through at (21/2,12,0), direction (1,-1,1), flat tip. Observed missing exact rational coefficient admission for traced cone/cylinder native window. Radial control is the same collar with point (9,12,0), direction (0,-1,0), matching the public success test. The control changes placement and direction and is not proof that every skew hole fails.

All derivations are implemented by [make_inputs.py](../../../evidence/phase1/P03/make_inputs.py), using only the baseline's public files and explicit selections above. Generated inputs retain source names, which can repeat across cases; use input filename and case key as identifiers.

## Replay and retained artifacts

From the isolated repository root, set `PROOF` to a new authorized scratch directory. The command examples below use a portable relative binary location; supply the preserved baseline binary identified above:

```sh
python3 evidence/phase1/P03/make_inputs.py
python3 .claude/skills/verify-catallax/scripts/control-catallax.py doctor --binary target/release/catallax-verify --json
python3 .claude/skills/verify-catallax/scripts/control-catallax.py run evidence/phase1/P03/inputs --source recipe --binary target/release/catallax-verify --out "$PROOF/recipe" --max-work 20000000 --workers 1 --cores 0 --fresh --json
```

For a new replay use a new external output directory and an available core reported by doctor. Each preserved `question.json` retains the original 20M allowance and is also replayable via `catallax-verify replay QUESTION --out NEW_DIR`. Replay/determinism comparison was not separately executed in this wave.

[raw](../../../evidence/phase1/P03/raw) contains **sanitized publication copies**, despite the retained directory name: doctor output, run stdout, configuration, ordered inputs, complete run, cases, and cleanup dry-run/actual output. Machine-specific locations have been replaced by `<PROOF>`, a relative binary location, or a portable description. These copies are not original raw bytes. Original evidence remains privately preserved in existing authorized scratch. No geometry execution was repeated for this correction. The archived files are under `raw/recipe/cases/KEY`; keys remain original execution identifiers and are not hashes of sanitized directories. Per-case artifact manifests and [SHA256SUMS.json](../../../evidence/phase1/P03/SHA256SUMS.json) were recomputed for the published bytes; previous hashes are not asserted to cover changed bytes. Binary and unchanged geometry-file hashes retain their original meanings. Original execution artifact checks, input/problem mapping and all 11 outcomes were checked before sanitization. Cleanup preserved evidence and no worker remained owned by the driver. See [publication notes](../../../evidence/phase1/P03/PUBLICATION.md) for the precise transformation.

## Limits and remaining requirements

This covers baseline CLI reproduction, not SolidDraft user flows, edit/save/reopen, display, native field, sampled comparison, requested mass properties, independent STEP import or a whole-result geometric qualification. No companion repository, private/held-out inputs, C++ oracle, credentials, hosted gates, CI or memory/sanitizer probes were used. No need for a memory-fault probe emerged from these completed typed outcomes. A routine attempt to list `a restricted toolchain location` while locating a toolchain returned filesystem permission denied; it was not retried or bypassed and no build relied on that path.

A fresh source build requires an available Rust toolchain. Separate Fillet/Wall/Bevel decision records, if Firstmate intended documents beyond these public sources, require an exact visible input location. Full application and independent consumer qualification remain outside this card's executed evidence. Agent-generated report and input extraction were checked against retained CLI artifacts; no algorithms, adapters, PRD, tolerances, allowances or saved/proven answers were edited.
