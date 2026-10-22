# Verification tooling checkpoint, 22 October 2026

This is the first target-part report, not a C1 completion claim. The [machine-readable report](2026-10-22-parts.json) was generated from commit `d98ac6a00552760399e2cd5f58b55a4a90179c0b`, after rebasing onto the owner's published roadmap, with binary SHA-256 `6dc829c48dc39381dd4557afb8b4ef78866142ae189c3f57894a9f91393b26ff`. Recipe hashes, complete adapter responses, work counts, operation logs and timings are in that report.

The per-request allowance was 1,000,000 units. The worker was pinned to CPU 2 on a shared machine whose initial 1/5/15-minute load averages were 20.18/21.16/20.35. These are diagnostic timings under load, not evidence of quiet-core speed or a slowdown. The B-rep adapter's time includes topology summary and STEP export when they run.

| Observation | Count |
| --- | ---: |
| Complete target recipes attempted in all three adapters | 13 |
| B-rep builds with STEP export | 3 |
| Independent field builds | 9 |
| OpenCASCADE imports passing validity and face-count checks | 3 |
| Proven sampled comparison disagreements | 0 |
| Internal errors or structurally invalid recipes | 0 |
| Completed target parts under all roadmap acceptance gates | 0 |

P0's equal-radius cylinder/sphere, line-tangent cylinders and near-miss cylinders export STEP that opens in OpenCASCADE 8.0.1.1 with valid topology and matching face counts. Their imported volumes are recorded. Kernel volume comparison is explicitly unmeasured; no volume-deviation claim follows.

The remaining P0 builds either exhaust the fixed allowance or report `intersect.traced-pair-required`. P1's B-rep build exhausts the allowance, while its independent field build completes. P2–P5 report `features.unsupported-operation` and `field.unsupported-operation` at real chamfer, fillet, shell or pattern steps. No unsupported part was omitted or replaced by an easier model. No comparison is reported as agreement; each is unresolved or exhausts the allowance.

The tooling now writes versioned per-operation logs and portable application questions. A 4×5×6 block recipe replayed through the exported question with the same manifest, 2,400 work units and byte-identical STEP. The permanent replay test also covers an unsupported feature, its stable refusal code, and the absence of later unexecuted steps. Entry-checkpoint work refusals retain their operation in the log.

The report runs on `main` through the separate `target parts` workflow and uploads its full artifact bundle. Its fixed dimensions are documented in [parts/README.md](../../parts/README.md); the schemas and coverage limits are in [VERIFY-INTERFACE.md](../VERIFY-INTERFACE.md). Application edits, persistent-name retention, kernel-volume comparison and whole-result tolerance checks remain separate acceptance gates. This checkpoint does not claim they ran.
