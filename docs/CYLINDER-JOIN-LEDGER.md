# Bounded cylinder join lane

Status: development; draft only. No readiness, replay-477, or timing qualification claimed.

## Frozen base and ownership

Fetched main/base: `71f552f5f3488601b8d238b118c9979a23f5a398` (2026-10-25).
The saved `work` branch was clean at `4b7ce8f0632feff5e327fbd1aeda164657dc640a` and is preserved.
Working branch: `cloud/offset-cylinder-joins`.
No AGENTS.md was found in the repository or its workspace ancestors. Read the repository verify-catallax skill, feature map, answer instructions and design/architecture contracts.

Active PR inspection before shared edits: #103 `759e83db11bac7297af173ccfda1dc235133be45`; #108 `ebfd11e3ce71d99a9e4b9bd8a30d2d729191e0e3`; #93 `4b18be87b0af640c9de18ece3e71354eb868a5b4`; #101 `baabe04998300e93fd28933999b584b59f9de394`. Connector metadata and patches inspected. PR103 remains unmerged; integration is required before final readiness.

Owned shared hunks, declared before editing:

- `catallax-features/src/lib.rs`: cylinder_join public module and private cylinder_join_build module registration next to cylinder_skin only.
- `catallax-field/src/lib.rs`: private cylinder_join module registration next to circular_distance only.
- `catallax-features/src/distance_build.rs::offset`: positive bounded join dispatch before circular_build recognition; no other offset behavior changes.
- `catallax-features/src/distance_build.rs::thicken`: bounded selected join face dispatch before cylindrical recognition.
- `catallax-field/src/distance_features.rs::distance_offset` and `distance_thicken`: equivalent bounded join dispatch at the same feature stage.
- New cylinder-join definition/construction/field modules, fixtures, tests and this ledger.

No ownership of triangulate.rs (especially PR103's per-triangle vmax hunk), display/P3, features.rs/convex_revolve, revolution, chamfer, Boolean machinery, number algebra or edge_blend. No elliptic-prism carrier, fitting, tolerance or work-budget changes. User-authorized 100 nm target and 1 um hard cap remain unchanged.

## Evidence and defects

Evidence outside the repository: `/workspace/cylinder-join-proof`. Public authored fixtures only; no private corpus or held-out input access.

An initial attempted baseline reproduction inadvertently used the preserved stale binary (`9e0200c`) while the fresh base build was still running. Its recipe/4 parse errors are invalid as geometry evidence, retained separately; repeat with the completed base executable before production edits.

Pending: fresh-base refusals; exact finite-domain construction and field witnesses; volume/point/display/STEP controls; negative domains; deterministic runs; focused then full final-head checks; PR103 landing integration; CI; draft publication.

## Development checkpoint (2026-10-25)

Fresh base executable `base-71f552f-verify` has SHA-256 `73878365b28c2e4a204fdddfa0b0eafb82f70a419dc3b7236e0557e81304f154`. Doctor confirms `catallax-request/1` and base identity. Commands `control-catallax.py run crates/catallax-features/tests/fixtures/cylinder-joins --source recipe|field --binary /workspace/cylinder-join-proof/base-71f552f-verify --out /workspace/cylinder-join-proof/base-71f552f-brep|field --json` reproduced four refusals each: B-rep `distance.cylinder-domain`, field `distance.unsupported-envelope`. Default verification allowance was 100M; this is not the 1M scoreboard or 20M application allowance.

The source witness checks every circular section, opposite incident uses, oriented patch support and complete rectangular outer segment before B-rep construction. Field recognition uses only complete primitive domains at the two-operand union. No arbitrary trimmed leaf is inferred from supports. A boss must attach perpendicularly and remain strictly contained; a shaft has two coaxial consecutive unequal-radius bands. Strict radial, axial, shoulder and plate clearances exclude new contacts. General trimmed/spline envelopes, noncoaxial/oblique joins and new sharp multi-face selections remain `distance.unsupported-envelope`.

Outward convex circular turns get exact torus sectors; concave turns intersect exact displaced plane/cylinder supports. The plate reuses rounded_polyhedron's edge cylinders and corner spheres. One selected smooth wall, circular disk/annular shoulder, or rectangular plate face with its circular hole gets its exact finite normal image and cylindrical rims. Existing supported stock/circular/conical and opt-in elliptic paths are not rewritten.

Focused `cargo test --release -p catallax-features --test cylinder_joins -- --skip field_volume`: 8 passed (1.95s), including exact volume enclosures, both-kernel point controls, negative clearances/noncoaxial/multi-face cases, selected wall/annular/plate trims, posed/scaled copies and closed display meshes at 0.005 mm chord control. This is a development dirty-tree result, not final-head qualification.

All four authored B-rep CLI fixtures answered and wrote STEP. `CATALLAX_CYLINDER_JOIN_EVIDENCE=/workspace/cylinder-join-proof python3 -B -m unittest discover -s tools/step-regressions -p test_cylinder_joins.py -v`: passed (0.653s). OpenCASCADE 7.8.1.1 imported all four without repair, with 29/8/4/4 faces, independent numeric volume readback and material controls. Import validity/numeric readback is separate from exact enclosure/deviation certification.

Independent analytic volumes at d=1/4 mm:

- Boss on 20x16x3 plate: `1174 + (1761/64) pi + (1/16) pi^2` mm3.
- Two-band shaft (R4 h3, R2 h5): `(5585/64) pi + (5/16) pi^2` mm3.
- R2 h5 outward wall normal sweep: `(85/16) pi` mm3.

Development defects retained in evidence: an initial plate assembly lacked a curve/surface incidence fact (fixed through the existing face constructor); the inherited field decomposition gave indeterminate internal seams (new normal pieces now overlap only in proved source interior); unreduced rational expressions prevented rational native charts (exact reduction in the new module only). Comparison/display reruns are pending.

Default field mass integration has an enabled failing regression. The initial boss request exhausted the unchanged 2M-cell limit (`field.mass-budget`, 51.97s); no mass accuracy or work/cell limit was enlarged. The later representation is being rechecked at the 100M verification allowance. A bounded exact mass hook in field `mass.rs::properties` plus retained join metadata in Plan is an additional ownership request, not an applied shared-file edit. Other mass/display limits remain visible.

One obsolete native-display development run was interrupted after its representation changed; its incomplete request is not a geometric refusal or latency result. Completed wall-thicken displays retained 0.01 mm native boundary certificates. All wall times here are cloud diagnostics, not pinned matched-pair speed qualification. Build time, operation work, independent runs and cache reuse remain separate.

Routine prerequisites: saved Rust paths supplied explicitly; missing m4 restored under the external proof directory, Debian trixie 1.4.19-8, SHA-256 `221b3e224708c1da2a1a88b6824380593579b2254ae7d914afa14cbd5d343334` verified against its package index. CLI GitHub GraphQL/REST reads returned Forbidden; the connector read PR metadata/patches successfully. No automatic approval rejection occurred.

PR103 remains open at `759e83db11bac7297af173ccfda1dc235133be45`. No merge, cherry-pick or duplication of its hunk has occurred. Landing integration and all final-head/CI evidence remain pending; this branch must stay draft.

## Frozen runtime qualification and handoff

Frozen implementation/runtime head: `dc434300dc0bb8f3379fb2d9461fab8bdeb577fb`, integrated on `1337f0ff3c44f1fb22454f36d9055689ed2f11f5` by normal merge. Immutable verification executable SHA-256: `889572b0524c1576bddcd3bc609001b5aef457d046f63fc1a0d445b7324b56a7`. Full commands, raw responses, manifests, STEP/OBJ, diagnostics and exact owned hunks remain under `/workspace/cylinder-join-proof`; `defect-ledger-dc434300.json` is the machine-readable handoff ledger.

- `cargo test --workspace --release --no-fail-fast -j 2 -- --test-threads 2`: completed at the frozen runtime head, **1,073 passed, 3 failed, 11 ignored**, including doctests. The three failed targets are CLI exchange and the two new field gates. A concurrent normal CLI build replaced the dev-only `test-deadline` binary during the initial run; its marker was absent, explaining the extra exchange assertion. Restoring the proper test-feature binary makes the focused deadline test pass in 0.10s, matching the independently built main-base test and CI. No source was changed. The original failed full run is retained and is not relabeled green. Its 40m11s wall includes compilation and execution under changing load, not latency qualification.
- [Exact runtime-head check CI](https://github.com/BohnBawerick/catallax/actions/runs/40236728498): **1,072 tests passed, two failed**. The only failures are `catallax-features::cylinder_joins::field_volume_enclosures_contain_the_same_exact_integrals` and `catallax-field::cylinder_joins::bounded_cylinder_join_native_display_has_complete_certified_outputs`. Both other test partitions, build (including native STEP and doctests), format and Clippy pass. [Target-parts CI](https://github.com/BohnBawerick/catallax/actions/runs/40236728511) passes. The legacy deadline test passes in CI.
- Fresh local runtime-head driver tests: **22 passed**. All native STEP regressions: **6 passed**, including all four new cases with expected 29/8/4/4 faces, independent volume readback and material controls. Format, Clippy with warnings denied, 15 edge-hunt tests, 21 release-policy tests and generated-page check pass.
- All four new recipes pass integrated two-kernel comparison with no findings. Local full-suite comparison took 68.97s; CI took 65.071s. These are shared-machine diagnostic observations, not timing qualification or whole-result distance proofs.
- Native field display at the 100M verification allowance: shaft offset answers with 3,603,117 work and 121,240 cells; boss/shaft wall thickens answer with 31,878/20,515 work. The boss offset refuses `field.display-budget` at 2,000,001 cells after 10,815,051 work. Rounded-plate straight-edge carriers dominate refinement. The default 0.01 mm display tolerance and all limits remain unchanged. A fresh one-worker CLI run observed 201.8s for the boss refusal and 57.3s for the shaft answer: raw misses against the unchanged 1000ms field click target, not an accepted replacement baseline.
- All four default-accuracy field mass requests refuse `field.mass-budget` at the unchanged 2M-cell limit. Their exact analytic oracle assertions remain enabled. The additional ownership request is a bounded dispatch in field `mass.rs::properties` and retained `Plan` metadata, with integration implemented in a new cylinder-join module. No such shared edit has been taken. Shared display classification/coverage changes would overlap PR108 and remain with that owner.
- Fresh independent one-, two- and four-worker CLI runs match all four responses in both sources, including deterministic work/decision counters and strict STEP/OBJ bytes. The expected boss-display refusal matches too. Four cache reuses per source match and are recorded separately from fresh execution. Zero work refuses and leaves no STEP.
- Against an independently built `1337f0f` executable, all **27 existing public distance fixtures** match in both evaluator entries (54 manifests/outcomes); the 25 B-rep answers and two earlier conical refusals are unchanged, as are the field's 15 answers, two conical refusals and ten work refusals at its separate 1M control allowance. The public **17-row 1M scoreboard** additionally preserves 51 manifests, complete work/operation logs and edit checks, with all seven STEP and five OBJ files byte-identical. Its earlier refusals remain refusals. This is public regression evidence, not private replay477 or complete target-part acceptance.

Current upstream at handoff inspection is `e75c2796bcdbc30ce03c2d8042e477375f539f3e` (PR99, PR110 and PR111 after the frozen runtime base). These upstream runtime mesh changes are **not** silently included in the frozen qualification. PR103 is still unmerged, now at `c2df33e5a806e9147f2dfa39acd3acdc36dc313e`; its per-triangle vmax hunk remains untouched. Further main/PR103 integration and fresh applicable qualification are readiness gates.

The final documentation update reads PR111's fragment instructions and adds only this lane's new interface fragment; there was no lane entry in the aggregate to remove. `docs/INTERFACE-CHANGES.md` is untouched. A before/after hash verifies the crate code patch is identical across this documentation update. Frozen runtime evidence is preserved under its exact head and is not relabeled as execution on the later documentation head. Final delivery identity and current CI are recorded in the PR description.

Status remains **DRAFT / NOT READY**: field mass, boss native display, PR103 integration, final qualification and parent replay477 remain open. No private proof-input contents or held-out inputs were read or used. No merge, readiness transition, force push, protections bypass or tolerance/budget change has occurred.
