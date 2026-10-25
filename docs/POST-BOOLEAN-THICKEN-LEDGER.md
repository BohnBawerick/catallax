# One post-Boolean polygonal face thicken

Base: `357be43e3ea5ba00df9fd449ea2afca09260ae24`. This lane owns only the two
`planar_face.rs` modules, the field module's `recovery.rs` and `sweep.rs` helpers,
and new dedicated tests, public fixtures and this additive documentation.
The unchanged dispatch calls are features `distance_build.rs:503` and field
`distance_features.rs:731`. No display, Boolean, shapes, numbers, shared
`features.rs`/`lib.rs` or held production file is changed.

The branch is `cloud/post-boolean-polygon-face-thicken`. Exact published head,
binary identity, commands, outcomes and pending final gates are recorded in the
PR description and `/workspace/post-boolean-thicken-proof`. Raw evidence stays
outside the repository. This is a functional lane report, not full target-part
acceptance or private qualification. final477 remains parent-routed.

## Reproduction and provenance

Before production edits, current main and the saved clean checkout agreed at
the base above. Open PR ownership and complete filename lists for
93/101/103/107/109/112 were inspected through the connected GitHub app. No
`AGENTS.md` was present at the filesystem/workspace/repository roots or below
the repository. The verification skill, feature map, deterministic-answer
entry, workflow checks and interface-fragment instructions were read.

A fresh release verifier at that exact base rejects all seven initial fixtures
with `distance.unsupported-envelope` in both evaluators. The B-rep rejects the
non-plain-prism body; field recovery rejects its non-leaf Boolean tape.
The preinstalled verifier was older (`9e0200c`) and rejected recipe/4; those
invalid initial attempts are retained separately as `stale-preinstalled-*`.
They are not geometric reproduction evidence.

All fixtures are new public synthetic definitions. The L is a 4 by 4 by 3 block
minus its 2 by 2 through-corner tool: selected area 12. The hole is a 6 by 6 by 3
block minus a 2 by 2 through tool: selected area 32. Unit-thickness front/back
sweeps therefore have exact volumes 12 and 32; the explicit front/back L union
has volume 24. Additional fixtures cover a reversed tool wall (volume 6), split
origin (volume 8), merged cap aliases and an exact rational rigid tilt.
No private corpus or held-out input was accessed.

## Construction and proof boundaries

The B-rep reads its current selected face loops, orientation and adjacent
physical face origins. The field independently collects planar source skins
from its active leaf definitions, classifies full-dimensional arrangement cells
through the Boolean tape and proves both outward germs and surviving origins.
It never calls the B-rep evaluator or uses native field display as recovery.

Both then construct the finite straight normal sweep. Source and offset caps
retain every alias; each rim retains the product of selected and adjacent
origin sets. Field convex material terms use exact analytic leaves. A term is
reduced to a stock box only when every bounding plane proves that exact box;
other terms keep all analytic halfspaces. Only physical boundary supports get
face candidates. Partition walls carry no invented origins. Primitive sampling
roles remain separate from public construction names.

The bounded field fallback receives a node, not recipe ownership snapshots.
It refuses geometry-folded creators with unrelated aliases, copied/derived
sources, curved source skins and multiple physical selections rather than
guessing their scoped provenance or normal fans. Lost/ambiguous selections
stay typed. Existing authored plain-cap behavior uses its original path.
The field registry retains historical candidates by design; physical removal
of the two-sided base is established by `face_query` and interior membership,
not by registry-list absence.

## Defect ledger

| Attempt | Finding | Disposition |
| --- | --- | --- |
| Setup | cargo absent from PATH; default rustup home read-only | Use existing `/workspace/.cargo` and `/workspace/.rustup`; no account/configuration change. |
| Frozen base | seven geometric refusals in each evaluator | Reproduced before production edits and preserved with binary identity. |
| First candidate | polyhedral Draft leaf refused field mass with `field.print-coverage` | Replace only the new helper's publication with certified analytic material terms. |
| First candidate | union name-set assertion counted a historical base candidate as physical | Test its actual trim membership at former base interiors and require inside material. |
| Analytic candidate | sampling `rim` names had no primitive chart (`field.sample-role`) | Preserve exact primitive definition roles separately in the owned helper. |
| Local Clippy setup | gmp build prerequisite `m4` absent | Repair only workspace-local prerequisites; final result recorded in PR/evidence. |
| Native STEP harness | OpenCASCADE 8 renamed `TopoDS.Solid_s` | Use the binding's native solid cast; nine topology/volume/material read-backs then pass without repair. |

Eleven focused geometry/oracle/edit/name/refusal tests and the nine-fixture
sampled-comparison test pass in the development snapshot. Nine independent
OpenCASCADE 8.0.1.1.0 STEP read-backs pass topology, numeric volume and material
controls without repair. The 29 retained-domain tests pass, including original
nonconvex and sharp thicken behavior; a separate plain-cap collinear-alias test
keeps its original naming choice. Final-head full checks, CLI preservation and
CI results must be read from their actual evidence, not inferred from this
checkpoint. Earlier failed attempts are not removed or relabelled as passing.

## Commands and limits

Rust commands use `CARGO_HOME=/workspace/.cargo RUSTUP_HOME=/workspace/.rustup`
and prepend `/workspace/.cargo/bin` to PATH; builds use `nice -n 10` and `-j 2`.
The exact base build is `cargo build --release -p catallax-cli --bins -j 2`.
Doctor and independent recipe/field runs use the repository's
`.claude/skills/verify-catallax/scripts/control-catallax.py`, the frozen baseline
executable and separate output directories under the external proof root.
Focused commands are `cargo test --release -p catallax-features -p
catallax-compare --test boolean_face_thicken -j 2` and the retained
`nonconvex_thicken`, `distance_features`, `sharp_thicken` integration targets.
Final commands and exit codes are saved verbatim in that proof root.

No work allowance, tolerance, approximation bound, cell cap, comparison option
or typed fail-closed policy changes. The 1M scoreboard, 20M application and
100M verification allowances remain distinct. Sampled comparison is not a
whole-result distance certificate; STEP read-back is not certified deviation.
No matched latency experiment or speedup claim is made.
