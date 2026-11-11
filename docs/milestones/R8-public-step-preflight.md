# R8 public STEP census: acquisition blocked

Date: 11 November 2026. Status: **preflight only; R8 remains incomplete**.
No real-part census, control import, or compatibility measurement ran.

The inspected current-main baseline is
`eee8fd123ce8a97e1115013d436f9db7e51d3043`. No census binary was built or
executed, so there is no measured binary identity. Unmerged M24 changes were
not included. The accompanying [provenance manifest](R8-public-step-provenance.json)
records the candidate archive and empty selection explicitly.

## Candidate and licence evidence

The [ABC authors' site](https://deep-geometry.github.io/abc-dataset/) describes
the April 2019 v00 STEP release, divided into chunks of 10,000 models.
[NYU's chunk 0000 record](https://archive.nyu.edu/handle/2451/44309) identifies
the 2019 archive and labels it MIT License. The
[author repository licence](https://github.com/deep-geometry/abc-dataset/blob/392789113158b874bfa630098f5bbab3ed38309e/LICENSE)
is MIT, copyright 2019 Deep Geometry Processing: research use is permitted;
copies or substantial portions must retain the copyright and permission notice,
and the software is provided without warranty. No CAD files are redistributed here.
Dataset attribution: Koch et al., *ABC: A Big CAD Model Dataset for Geometric
Deep Learning*, CVPR 2019.

The authors and NYU also say original CAD copyright remains with creators and
refer to Onshape terms. The linked
[Onshape terms, effective 15 July 2020](https://www.onshape.com/en/legal/terms-of-use),
section 1.g.ii, grant broad use of qualifying Free Plan public documents but
retain exceptions for some documents with their own licence tabs. The current
terms postdate the dataset release. The dataset-level MIT label is useful
provenance, not verification of every source document's rights. Per-model
metadata and applicable exceptions must be resolved before accepting files
into an eligible sample. Public visibility alone is insufficient.

Author-published size metadata gives `abc_0000_step_v00.7z` as **1,594,129,754
bytes** (about 1.59 GB decimal); published MD5 is
`695388be7a278c7798c8c8ae239772ac`. This is one chunk, not the full corpus.
The manifest pins the source-list and size/checksum metadata versions. These
are expected archive properties, not a verified downloaded archive checksum.

An alternative was investigated but not downloaded:
[Fusion360 Gallery](https://github.com/AutodeskAILab/Fusion360GalleryDataset)
advertises an extended segmentation STEP release with 42,912 files (483 MB).
Its [November 2021 custom licence](https://github.com/AutodeskAILab/Fusion360GalleryDataset/blob/master/LICENSE.md)
restricts use to non-commercial research, restricts redistribution, requires
notices/attribution, and includes employer-authority conditions. Applicability
of those conditions to this run was not established; it was not substituted
for ABC or described as MIT/CC-BY.

## Observed acquisition blocker

A bounded `curl --head --location --max-time 30` probe to `archive.nyu.edu`
failed at the executor's proxy CONNECT stage:

```text
HTTP/1.1 403 Forbidden
curl: (56) CONNECT tunnel failed, response 403
```

The response date was `Wed, 11 Nov 2026 04:30:44 GMT`. The probe used an
unverified `/bitstream/2451/44309/31/abc_0000_step_v00.7z` path. Subsequently
the authoritative record established `/3/` and the author source list gave
`https://archive.nyu.edu/rest/bitstreams/88598/retrieve`. Because CONNECT was
denied before contacting the origin, the result establishes an executor
host-access blocker, **not** an unavailable archive, invalid endpoint, or
licence refusal. No alternate download route or repeat denied probe was used.
Only small public GitHub licence/size/checksum metadata was fetched.

Completion needs permitted executor access to the archive host, followed by
the remaining per-source licence and archive-layout checks. No account or
credential change was made or requested by the tool. This report does not
claim that access alone completes the licence review.

## Coverage and denominators

| Quantity | Observed result |
|---|---:|
| Intended distinct eligible sample | 3,000 |
| Selected / downloaded / evaluated STEP parts | 0 / 0 / 0 |
| Input hashes and duplicate groups | None collected |
| Referenced faces and per-kind frequencies | Not measured |
| Parts whose referenced face kinds are all buildable | Unknown; denominator 0 |
| Full import admission / typed refusals / parser failures | Not observed |
| Per-part timeout or resource exclusions | Not observed |
| Acquisition failures | 1 host CONNECT denial |
| Fillet/chamfer counts | Not measured |

Compatibility ratios are null, not zero percent. The 3,000 intended parts
are not 3,000 failures or resource exclusions: they were never selected.
No successful parse, first refusal, or surface tag has been counted as a
complete solid import. There are no measured refusal diagnostics to tabulate.

## Resume protocol and measurement boundaries

1. Check the published archive size/layout before transferring data. Prefer
   an indexed subset if the archive supports it; do not assume 7z permits
   independent remote extraction. If a whole chunk is necessary, assess the
   single 1.59 GB compressed chunk and decompression/storage cost first; do
   not fetch the full corpus. Retain source licence metadata for eligibility.
2. Sort eligible model identities bytewise within chunk 0000 and select the
   first 3,000 distinct content SHA-256 hashes, recording each selected entry,
   duplicate, exclusion and shortfall. Deduplicate exact bytes only; geometric
   duplicates remain a disclosed limitation. This deterministic chunk sample
   would not be an unbiased estimate of the million-model corpus.
3. Build a Rust-only census against an explicitly recorded current-main SHA,
   retain compiler/build options and executable SHA-256, and keep its output
   separate from production decisions. No census tool has been implemented
   or validated in this blocked attempt; these are requirements, not runnable
   commands or completed work.
4. Traverse referenced solid/shell/face/trim/surface entities rather than
   count every surface token. Deduplicate shared entity references for source
   counts; separately record instances and transformations, unresolved
   references, unknown representations and incomplete coverage. Unsupported
   rational/analytic variants must remain distinguishable from supported
   carriers. A first reader refusal must not truncate the face census silently.
5. At the inspected baseline, `read_step` recognizes more than certified
   construction accepts. The [import contract](../pieces/08-step-import-repair.md)
   distinguishes analytic plane/cylinder/cone/sphere/torus carriers and their
   restricted trims/topology from general spline recognition: spline solid
   faces/edges are refused by `build_file`/`build_solid`. An analytic tag alone
   does not establish admissible trims, cone bounds, manifold topology or
   representation completeness. Use a fresh model per evaluation and the
   actual recipe import boundary; partial builder mutation is not admission.
   Record full diagnostics and repair requirements without auto-accepting
   repairs. Keep original work/resource controls; wall-time termination is
   a separate observation, never permission to increase a proof allowance.
6. Report both per-face and per-unique-part frequencies, all-kind coverage
   only for completely traversed parts, and separate admitted/refused/read
   failure/resource/unknown populations. Publish every denominator. Surface
   curvature does not establish fillet/chamfer feature provenance; report
   those unknown unless independently justified metadata is available.
7. Validate arithmetic, replay and complete traversal on a small public
   control set before the census, then retain a source/hash selection manifest
   and reproducible commands/results. Clean disposable downloaded copies.

This preflight changes no production code, answers, allowances, tolerances,
names or acceptance. It provides no M21/M24 compatibility conclusion or
kernel-correctness evidence. No Rust/C++ tests, oracle, saved-answer corpus,
private inputs or hosted CI reruns were used. Documentation/manifest checks
are the only validation appropriate to this partial deliverable.
