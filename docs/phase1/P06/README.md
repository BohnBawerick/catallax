# P06 — Checkpoint 1 comparison input pin

This is a partial input-definition report, not a Checkpoint 1 qualification.
Baseline: `BohnBawerick/catallax` main
`f4db345f2e7c51835d1ebeba9f6446d117144117`, verified against origin before
creating isolated branch `phase1/p06-checkpoint-inputs-20261114`.

[inputs.json](inputs.json) is the single ordered list: P1–P10 first, then ten
public census controls in the existing M5 answering order. It records 29 source
file SHA-256 hashes, exact committed edits, original operation IDs/order, source
licence scope and missing definitions. [pin_inputs.py](pin_inputs.py) reproduces
it from this baseline using Python standard-library text parsing and hashing.
Run from the checkout with `python3 docs/phase1/P06/pin_inputs.py`.
The generator does not import or execute repository kernel/scoreboard modules.

## Pinned edit intent

Each chain starts with the complete original target, followed by the complete
edited target with all IDs, origins, placement and nonedited parameters retained.
The same scalar change is intended in the reference kernel. Existing committed
reference commands are recorded, together with the reference release's published
archive/binary hashes and Apache-2.0 licence. Neither installed binaries nor
release bytes were inspected or downloaded. Reference definitions and licence
metadata are source evidence, not executed equivalence or current upstream checks.

| Target | Committed scalar change | Definition status |
| --- | --- | --- |
| P1 bracket | `boss.radius`: 8 → 9 | Full original and one edit pinned |
| P2 collar | `collar.distance`: 1 → 0.75 | Full original and one edit pinned |
| P3 filleted boss | `junction.radius`: 2 → 2.5 | Full original and one edit pinned |
| P4 shelled tray | `tray.thickness`: 2 → 1.5 | Full original and one edit pinned |
| P5 flanged housing | `filleted.radius`: 2 → 2.5 | Reference-case edit pinned; absent from `parts/edit-checks.json` |
| P6 lattice bracket | Undefined | Full recipe exists, including lattice, cuts and unions; comparison edit and reference lattice definition absent |
| P7 field-driven bracket with gyroid infill | Undefined | Missing fixed gyroid field/region, grading parameters, placements, operations and edit |
| P8 helical spring and threaded boss | Undefined | Missing helix/thread profiles, pitch, sweep construction, joins and edit |
| P9 lofted duct | Undefined | Missing sections, placements, loft correspondence/closure, wall construction and edit |
| P10 topology-optimised bracket | Undefined | Missing domain, constraints, loads/objective, optimisation/field construction and edit |

P7–P10 detail above names required decision inputs, not an invented recipe.
`docs/TARGET-PARTS.md` says later recipes are written when their turn comes;
P6 now exists in `parts/` despite that document's older grouping of the next five.
No committed continuation beyond one scalar edit was found for P1–P5. A longer
Checkpoint 1 edit sequence requires an owner definition; this report does not
promote five single edits into an established multi-edit benchmark.

Identical edit intent does not establish identical operation graphs. Existing
`targets.py` uses boxes for P1 extrusions and a capsule cut for the overlapping
slot chain, cylinders for P2/P5 revolved stock, and explicit P5 hole copies for
the pattern. The SolidCraft cases carry their own authored commands. Their
published intent is complete geometry, but equal geometry and semantic naming
have not been established here. An operation-by-operation comparison needs a
reviewed common mapping, particularly sketch/extrude/revolve, pattern, lattice,
field/gyroid, sweep/thread, loft and optimisation operations. Unsupported later
operations must remain unavailable, without replacement geometry.

## Census inputs and provenance

The ten identities are ABC `00000131`, `00000173`, `00000190`, `00000405`,
`00000511`, `00000605`, `00000632`, `00000633`, `00000695`, `00000821`.
No identity is missing from the existing M5 ten-answering selection. Candidate
`00000549` remains explicitly excluded because the historical complete request
refused at its fixed allowance; the public scoreboard retains all eleven
candidates. This selection is a capability control, not a representative sample.
It is not asserted to be an independently specified Checkpoint 1 census selection.

Source: existing R8 public census, ABC chunk 0000 STEP v00 (2019), seed
`20261111`. Published archive SHA-256:
`52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe`.
The manifest pins archive size/MD5 too, public census rows, and the historical M5
exported STEP hashes. **Exported STEP hashes are not source geometry hashes.**
Individual raw source hashes, runtime wrappers and a specific post-import edit
for every census part are unavailable in the inspected public sources. No raw
corpus, private evidence or archive-member path was accessed. Equal input bytes
and identical edit semantics in a future reference run remain requirements.

The repository's own documentation/recipes use its MIT licence. `NOTICE` retains
separate Apache-2.0 notices for adapted implementation portions; those portions
are not copied here. R8 reports the ABC dataset MIT grant (Koch et al., *ABC: A
Big CAD Model Dataset for Geometric Deep Learning*, CVPR 2019), with individual
source exceptions unresolved. This is not per-model licence clearance. Only
existing public identifiers, metadata and hashes are republished.

## Executed checks and limits

Executed: fetched origin/main and matched its full SHA with the remote head;
created an isolated worktree; parsed public JSON/CSV; checked all five reference
case recipe hashes; resolved each scalar path in the full recipe; generated
canonical edited-recipe hashes; checked 20 unique ordered identities; regenerated
the manifest byte-identically; and checked whitespace/diff scope. These are static
checks, not geometry tests. Source file hashes apply to exact bytes at baseline.

No kernel, reference CAD tool, C++ oracle, sanitizer, memory-fault probe, model
preparation/download, hosted gate or CI launch was executed. No algorithms,
Result/build adapters, PRD, tolerances, allowances or saved answers changed.
The original worktree and branches were preserved.

Missing for a fully executable Checkpoint 1 comparison: an authoritative longer
edit-chain definition if required; P6's edit/reference definition; P7–P10 full
recipes and reference operations; census per-source hashes/authorized public
inputs and ten exact edits/reference definitions; and reviewed operation mapping
where existing reference construction differs. Execution, correctness, naming,
application behavior, proof completion and speed qualification remain unmeasured.
