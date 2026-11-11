## 2026-11-11: a versioned seam-curve contract for STEP output

Pieces: 7, judge interface, verification exchange, CLI, STEP consumers

Change: STEP output gains an opt-in contract under which every curved face
carries its seam curves and pcurves, or the request refuses. ND5 still governs
the legacy exchange, together with the dated exceptions that already write
seams for complete bands, apex cones, loopless tori and regular partial charts
([band seams](2026-10-26-periodic-analytic-band-step-seams.md),
[loopless faces](2026-10-26-loopless-step-natural-boundaries.md),
[apex cones](2026-10-28-step-apex-cone-bands.md),
[partial charts](2026-10-26-partial-periodic-step-trims.md)). The new rule,
ND12, is reached only through a new version:

- `catallax_step::StepOptions` gains `seams: Seams`. `Seams::Legacy` is the
  default and writes exactly the bytes it wrote before. `Seams::Required`
  writes the same bytes, then checks the finished file before any byte leaves.
  `catallax_step::check_seams(&[u8])` applies the same check to a file written
  by any other writer.
- `exact-brep-question/3` has the schema of `exact-brep-question/2`. It is
  answered with `exact-brep-manifest/3`, which has the fields of manifest/2;
  its `step` object adds `"seams": "required"`. `exact-brep about` lists both
  question formats and adds `seamed_manifest_format`.
- The verification exchange accepts `/3` questions on the `question` source and
  a `catallax-step/1` wrapper on the `recipe` source: exactly `format` and
  `recipe`, where `recipe` is any problem the recipe source reads. Its `step`
  object adds `"seams": "required"`. Certified fitted results are checked after
  their own writer runs. `catallax-verify about` adds `seamed_step`. A wrapper
  with other fields, or sent to another source, is `invalid-question` with
  `cli.step-wrapper`.

ND12: every face whose surface is not a `PLANE` places each edge on itself
through a `PCURVE` of a `SURFACE_CURVE` or `SEAM_CURVE`. An edge that the face
uses twice is a `SEAM_CURVE` with a pcurve on that face for each use. A
`VERTEX_LOOP` bounds only a `SPHERICAL_SURFACE`, where the writer puts it on the
reader's own pole. Otherwise the export refuses `step.seam-unsupported`, writes
no file, and names the count of failing faces by surface type and the first
failing face entity. The check reads the written entities, so it covers every
encoder. It is structural; the pcurves' certified distances are unchanged and
are part of `max_deviation` as before.

Why: [PRD](../PRD.md) and PLAN M24 ask for seam curves in STEP for interchange
with other CAD tools, through a versioned change rather than by silently
widening ND5. Under ND5 a reader cannot tell whether a file's curved faces carry
their seams. Under ND12 it can, and a file that would need the reader to build
seams is a typed refusal instead.

Measured coverage, raw evidence retained privately: the 237 public recipe and
question fixtures in this repository give 195 STEP files. On every one, the
ND12 check agrees with an independent raw OpenCASCADE read (the pinned
`cadquery-ocp==8.0.1.1.0`, with no ShapeFix, SameParameter or tolerance
limiting, followed by `BRepCheck_Analyzer`): all 176 files that meet ND12 are
valid, and all 19 that do not are invalid. All 19 come from recipes and refuse
`step.seam-unsupported` under the new contract: eight cones trimmed by planes
or limits, two oblique elliptic extrusions (`SURFACE_OF_LINEAR_EXTRUSION`), and
nine offset or distance-limited parts mixing cone, cylinder, sphere and torus
faces, where one face outside the chart writers keeps every face of the file on
the legacy encoding. Wrapped in `catallax-step/1`, the other 207 of the 226
recipes give the same manifests apart from the added `seams` field, and every
written file has the same bytes as without the wrapper.

Seamed files import again. The importer reads `SEAM_CURVE` through its 3D
curve and drops each edge used once each way in one bound. A cylinder,
frustum, pointed cone, sphere, torus, through hole and blind drill tip written
under ND12 import with the same face kinds and an overlapping certified volume,
and the imported body writes an ND12 file again. Tilted copies can refuse for
causes the R8 census already lists, such as a cone axis of irrational length
after decimal rounding; those are not seam refusals.

Migration: struct literals of `StepOptions` add `seams: Seams::Legacy` or use
`..StepOptions::default()`. Consumers that need seam curves send
`exact-brep-question/3` or `catallax-step/1` and handle `step.seam-unsupported`.
Exhaustive readers of `exact-brep about` accept the added key. `Manifest` in
`catallax-cli` gains a `format` field that is not serialized as a member; set
it with `Question::manifest_format`. No legacy byte, manifest, saved replay
answer, tolerance or work allowance changes.

Known gaps: the 19 refusal shapes above need chart writers for general cone
trims, elliptic extrusions and mixed scenes; that is coverage work, and until
it lands those files refuse under ND12 rather than leaving seams to the reader.
