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
