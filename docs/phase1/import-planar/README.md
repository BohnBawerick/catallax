# Planar numerical source admission

The [numerical import contract](../../interface-changes/2026-11-15-numerical-step-import.md)
describes this kernel admission step toward ROADMAP 1.2.
It adds a separate numerical result after a recorded planar source repair and
makes the existing certified consumer entry atomic across all source solids.

The two fixed census identities `00000022` and `00000134` both refused strict
import at baseline `e9887bad203e351cd138b3eb8003e4202400cc05`.
At producer `67eabb12e18de9bc0463829c7721d1fffee68756`, each returns a numerical
body with six faces, twelve edges and eight vertices, under the stated 0.01 mm
local tolerance. The actual retained source-to-result upper bounds are at most
0.000000000001 mm. Zero or small residuals do not establish an exact result.

[The machine-readable controls](admission-controls.json) pin both source hashes,
the actual binary hash, charged work, check counts and outcomes. Original sources
come from the same checksum-pinned ABC chunk and seed `20261111` as the 3,000-file
census. Geometry and full source/response records remain private under the
[existing source and licence scope](../../milestones/R8-public-step-census.md#source-and-licence-scope).

Both complete import/edit/export workflows still refuse at the existing N03
writer's oblique line/frame limit. These are two successful source admissions,
not complete workflow qualifications or a new census import share. The separate
writer extension and actual edited STEP round trips remain necessary.

All 270 import all-target tests and the complete `ci/run fmt` gate passed on the
work PC at the producer revision. Nine focused numerical controls pass locally,
including source limits, crossings, multi-solid rollback, a rigidly rotated
oblique box and an actual pose rebuild with approximate disk STEP reread.
The unchanged N03 axis-aligned family supplies the successful export fixture;
the two real oblique files are retained as later export refusals in this record.
CI-selected import/reverse-consumer all-target Clippy passed before publication.
