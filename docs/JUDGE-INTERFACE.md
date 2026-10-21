# Kernel side of the judge exchange

This is the public exchange contract for the B-rep half of Catallax and an independent tester built without seeing the kernel's code. The kernel accepts a question file, answers it, and writes a STEP AP242 file and a JSON manifest. The exchange consists only of those files and the program's exit code.
This file specifies the kernel's input and output behaviour, not the tester's implementation or checks.

The program is `exact-brep`, built by the `catallax-cli` crate. Format versions are part of the contract: a change to either file format bumps its version and is logged in `INTERFACE-CHANGES.md`.
The current versions are `exact-brep-question/2` and `exact-brep-manifest/2`. Version 1 of each was replaced before any build wrote or read it.

`exact-brep`, `exact-brep-question/2`, `exact-brep-manifest/2` and the manifest's `kernel.name` value `exact-brep` are this contract's legacy names. The `catallax` binary is built from the same main and accepts the same commands. Rename these protocol identities only through a versioned contract change, never as a package rename.

ND1 to ND11 below record the existing exchange decisions.

## Units and numbers

- Lengths are millimetres. Angles do not appear in questions.
- Every number in a question is a JSON **string** in the exact rational grammar of piece 1, which has two forms and no others (ND2): a decimal with an optional exponent, such as `"1.1"`, `"-2.5e-3"` or `"007"`, or a fraction `p/q`, such as `"1/3"` or `"-6/8"`. The kernel parses it exactly (`"1.1"` is 11/10). JSON numbers are rejected wherever a geometric value is expected, so no value passes through a double.
- Every number the kernel writes into the manifest (bounds, mark bounds, `max_deviation`) is a decimal string, never a fraction. It uses the notation of piece 1's `to_decimal_bounds`: plain digits when the decimal exponent is from -6 to 20, otherwise one digit, an optional fraction part and an exponent with its sign (`"1e+30"`, `"2.5e-7"`). Both notations parse back exactly by the question grammar.
- Bounds in the manifest are rounded outward, so the true value lies between `lo` and `hi`.

## The question file

Format `exact-brep-question/2`. A JSON object with these keys; unknown keys make the question invalid.

| Key | Required | Value |
| --- | --- | --- |
| `format` | yes | `"exact-brep-question/2"` |
| `id` | yes | Question id, 1 to 200 characters from `A-Z a-z 0-9 . _ -`. Copied into the manifest. |
| `operands` | yes | Object mapping operand names to operand objects. At least one. |
| `steps` | yes | Array of step objects, possibly empty, applied in order. |
| `result` | yes | Name of the operand or step whose body is the answer. |
| `measure` | no | Object; default: measure everything at the default accuracy. |
| `export` | no | Object; default: write STEP with default options. |
| `limits` | no | Object; default: `{"deadline_ms": 50000}`. |

Names (operand names and step ids) are 1 to 64 characters from `a-z 0-9 _ -`, unique across operands and steps.

The `id` is an opaque label (ND8). The kernel copies it into the manifest and uses it as the STEP product name, and reads nothing else from it. A question carries no pose, family or cell label, and the kernel never guesses one: it answers from the geometry alone. Operand names and step ids are opaque in the same way.

### Operands

```json
{ "shape": "cylinder", "radius": "5", "height": "20",
  "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["10", "0", "-2.5"] } }
```

| `shape` | Parameters (all positive unless stated) | Local frame |
| --- | --- | --- |
| `block` | `size`: `[sx, sy, sz]` | the box `[0, sx] x [0, sy] x [0, sz]` |
| `cylinder` | `radius`, `height` | axis along `z` from `z = 0` to `z = height` |
| `sphere` | `radius` | centred at the origin |
| `cone` | `base_radius` (> 0), `top_radius` (>= 0, not equal to `base_radius`), `height` | axis along `z`; radius `base_radius` at `z = 0`, `top_radius` at `z = height` |
| `torus` | `major`, `minor` (`0 < minor < major`) | centred at the origin about the `z` axis |

`pose` is optional (identity by default). `rotation` is a quaternion `[w, x, y, z]` of exact numbers, not all zero and not necessarily of unit length (ND1). Integer quaternions are the usual form; any rational quaternion is accepted, and multiplying all four parts by the same non-zero number gives the same rotation. With `n = w^2 + x^2 + y^2 + z^2` it rotates by the matrix

```
        | w^2 + x^2 - y^2 - z^2   2(xy - wz)              2(xz + wy)            |
(1/n) * | 2(xy + wz)              w^2 - x^2 + y^2 - z^2   2(yz - wx)            |
        | 2(xz - wy)              2(yz + wx)              w^2 - x^2 - y^2 + z^2 |
```

acting on column vectors. `translation` `[tx, ty, tz]` is added after the rotation: a local point `p` lands at `R p + t`. Rational quaternions reach every rotation as closely as wanted while keeping the posed shapes exact; a question that needs a particular angle uses a rational quaternion close to it.
A rotation given as a 3 by 3 matrix is not accepted. Every rotation matrix with rational entries is the matrix of some integer quaternion, so the quaternion form reaches every exact pose, and the kernel never has to check that a matrix is orthogonal.

Cones and tori are accepted by the format now; until their pieces cover them, a question using them is answered with a refusal (`intersect.unsupported-pair` or similar), not rejected as invalid.

### Seams and poles

Every surface of an operand is parametrised in that operand's local frame, with reference direction `+x` and axis `+z` (piece 2's `ref_dir` and `axis`), and the pose carries the parametrisation with the shape (ND11).

- The seam of each periodic surface, where its angle parameter `u` is 0, lies on the half-plane `y = 0, x > 0` of the local frame. For a cylinder or cone it is the line of the surface in that half-plane; for a sphere, the half great circle from pole to pole through `(radius, 0, 0)`; for a torus, the tube circle centred at `(major, 0, 0)`.
- A sphere's poles lie on the local `z` axis, at `(0, 0, -radius)` and `(0, 0, radius)`.
- A torus has a second seam, where its tube parameter `v` is 0: the outer equator, the circle `z = 0` of radius `major + minor`.
- A boolean never replaces an operand's surface, so every face of an answer keeps the seam of the operand it came from.
- The STEP file places each analytic surface with the same reference direction and axis, so a reader sees the seam in the same place. The vertex piece 7 adds to a ring edge sits at its circle's parameter origin; for the rims of a cylinder or cone that point lies on the seam.

### Steps

```json
{ "id": "s1", "op": "subtract", "args": ["a", "b"] }
```

`op` is `unite`, `subtract` or `intersect`. `args` names two operands or earlier steps. Chains are written as several steps; `result` names the last.
The two names may be the same (`["a", "a"]`). With `X` the body they name, `X ∪ X` and `X ∩ X` are `X`, and `X − X` is empty, as the regularised set operations give. The step makes no decision, and its result carries `X`'s marks once.

### Measure

```json
{ "properties": ["volume", "area", "centroid"], "relative": "1e-12", "absolute": "1e-30" }
```

`properties` lists any of `volume`, `area`, `centroid` (default all three). `relative` and `absolute` set the requested width of each enclosure: at most `max(relative * |value|, absolute)`.

### Export

```json
{ "step": true, "curve_tolerance": "1e-7", "significant_digits": 17 }
```

`step` (default `true`) asks for the STEP file. `curve_tolerance` (millimetres, default `"1e-7"`) bounds the distance between the STEP file and the exact body (see "The STEP file"). `significant_digits` is an integer from 6 to 40 (default 17; a JSON number is allowed here, since it is not geometry).

### Limits

```json
{ "deadline_ms": 50000 }
```

Wall-clock limit for the whole answer, counted from the start of the program. Work in progress is abandoned when it passes, and the program writes a refusal manifest with code `numbers.deadline` and exits 0 within 1 s of the deadline. The default of 50 s leaves 10 s inside the 60 s outer limit below for that.

### A complete example

```json
{
  "format": "exact-brep-question/2",
  "id": "q-7f3a9c01",
  "operands": {
    "a": { "shape": "block", "size": ["10", "10", "10"] },
    "b": { "shape": "cylinder", "radius": "3", "height": "12",
           "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["5", "5", "-1"] } }
  },
  "steps": [ { "id": "s1", "op": "subtract", "args": ["a", "b"] } ],
  "result": "s1",
  "measure": { "properties": ["volume", "area", "centroid"], "relative": "1e-15", "absolute": "1e-30" }
}
```

The true volume of this answer is `1000 - 90 pi`.

## The program

```
exact-brep answer <question.json> --out <dir>
exact-brep check <question.json>
exact-brep about
```

- `answer` reads the question, answers it and writes `<dir>/manifest.json` and, when a STEP file is requested and the answer is a solid, `<dir>/answer.step`. It creates `<dir>` if needed and overwrites those two files. When it writes no STEP file it removes any `answer.step` already in `<dir>`. It runs single-threaded.
- `check` validates a question without answering it. It prints nothing and exits 0 when the question is valid; otherwise it prints the reason on standard error and exits 65.
- `about` prints one JSON object on standard output: `name`, `version`, `commit`, `question_formats` (array), `manifest_format`, `shapes` (array), `operations` (array) and `properties` (array).

Exit codes of `answer`:

| Code | Meaning | Manifest |
| --- | --- | --- |
| 0 | Answered, refused or invalid question | written |
| 64 | Command-line usage error | not written |
| 70 | Internal error (a panic or a broken invariant inside the kernel) | written with status `internal-error` when possible |
| 74 | The output directory or files could not be written | not written |

A refusal is a normal outcome: it exits 0 with status `refused`, like an answer or an invalid question. Every nonzero exit (64, 70, 74 or any other), a signal, a run past the outer time limit, or a missing manifest is a crash, even when a manifest with status `internal-error` was written. The kernel never decides by itself that a crash is acceptable.

### Running environment

These hold for every command (ND9).

- The program reads only the question file and writes only the files named above. Apart from its own executable it needs nothing from the machine: no network, no environment variables, no configuration or home directory, no temporary files, no locale. It runs in an empty sandbox.
- The judge may stop the process after 60 s of wall-clock time or when it uses more than 4 GiB of memory. The kernel's own deadline (default 50 s) refuses before the time limit. The kernel does not measure its memory, so a process stopped for memory is a crash.
- With the same build, the same question gives a byte-identical manifest, apart from `stats`, and a byte-identical STEP file, on every run and every machine. The one exception is an answer cut by the deadline: a run refused with `numbers.deadline` may differ from a run that finished.

## The manifest

Format `exact-brep-manifest/2`. One JSON object, UTF-8, keys in the order below.

| Key | Value |
| --- | --- |
| `format` | `"exact-brep-manifest/2"` |
| `question_id` | The question's `id`, or `null` if the question could not be read that far |
| `kernel` | `{ "name": "exact-brep", "version": "...", "commit": "..." }` |
| `status` | `"answered"`, `"refused"`, `"invalid-question"` or `"internal-error"` |
| `exact` | `true` when `status` is `answered` and `marks` is empty; otherwise `false` |
| `result` | When answered: `{ "kind": "solid" or "empty", "lumps", "shells", "faces", "edges", "ring_edges", "vertices", "faces_by_kind": {...}, "edges_by_kind": {...} }`; otherwise `null` |
| `properties` | When answered: `{ "volume": B, "area": B, "centroid": [B, B, B] or null }` with only the requested keys, each `B` being `{ "lo": "...", "hi": "..." }`; otherwise `null` |
| `step` | When a STEP file was written: `{ "file": "answer.step", "max_deviation": "...", "approximated_curves": n, "added_vertices": n, "added_edges": n }`; otherwise `null` |
| `marks` | Array of every mark on the answer, in order: `{ "code", "kind": "distance" or "angle" or "value", "bound": "...", "step": step id or null, "detail" }` |
| `refusal` | When refused, invalid or internal error: `{ "code", "step": step id or null, "stage": "intersect", "imprint", "classify", "keep", "sew" or null, "detail" }`; otherwise `null` |
| `stats` | `{ "wall_ms", "decisions", "filter", "refined", "exact", "coincident_exact", "coincident_separation", "coincident_fact", "not_exact", "refused", "max_bits", "work", "per_step": [ { "id", "wall_ms", "decisions", "not_exact" } ] }` |

Meaning of the fields:

- `faces_by_kind` keys are `plane`, `cylinder`, `cone`, `sphere`, `torus`; `edges_by_kind` keys are `line`, `circle`, `ellipse`, `intersection`. Kinds with a count of zero are left out. An edge counts by the shape of its curve: an edge where two surfaces meet in a line, circle or ellipse counts as that kind even when its definition is a meeting of two surfaces, and `intersection` counts only branches with no closed form.
- An empty answer is an answer (ND3): `status` is `answered`, `result.kind` is `empty` and every count is 0, `volume` and `area` are exactly `"0"`, `centroid` is `null`, `step` is `null` and no STEP file is written. `exact` and `marks` work as for any other answer. A refusal is never reported as empty, and an empty answer is never reported as refused.
- Every `lo` and `hi` is a proven bound for the body the kernel built. When `exact` is `false`, that body may differ from the true answer; each mark says where the kernel treated an unproven case as coincident and proves how far the true case can be from it (`bound`, in millimetres for `distance`, radians for `angle`, and in the decided quantity's own units for `value`).
- `properties` are computed from the exact definitions of the built body, never from the rounded STEP copy (ND10). Marks whose code starts with `measure.` or `step.` mean an enclosure is wider, or the STEP copy further away, than asked; they never change the body. So when every mark is of those two kinds, the built body is the true answer and each interval contains the true value of its property. When any other mark is present, the intervals are still proven for the built body, but the kernel does not claim they contain the true answer's values.
- `max_deviation` (millimetres) is a proven upper bound on the distance between the STEP file's geometry and the exact body.
- `refusal.code` is a stable `piece.reason` code from the piece descriptions, or one of `cli.invalid-question`, `cli.internal-error`. `detail` fields are free text and not part of the contract.
- `stats` is measurement only. It is the one key whose content may differ between runs.

### The answer's topology

The counts in `result` describe the answer body's own topology. Entities that piece 7 adds only because STEP needs them are reported in `step.added_vertices` and `step.added_edges` and are not in `result`.

Faces and edges are maximal (ND6):

- no edge separates two faces that lie on the same surface with the same sense;
- no edge has both its coedges in one face, and a vertex used only by such an edge goes with it;
- no vertex joins exactly two edges that lie on the same curve and bound the same two faces;
- an edge that covers a whole closed curve has no vertex unless another edge uses that vertex or the vertex is a singular point of the curve (the node of two meeting branches), so it is a ring edge.

Two faces on the same surface that share no edge stay separate faces. The kernel decides "same surface" and "same curve" with certified decisions; when one ends not exact, the boolean is refused with `boolean.merge-ambiguous`, because the merge would need facts that no proof gives. So an answered topology is always maximal.

Where the answer's material touches itself, the two sides of the contact share no face, edge or vertex (ND4). Each side keeps its own entities, which may lie on the same curve or point.

- Pieces of material that touch each other only along curves or at points are separate lumps. Two blocks united along a common edge give 2 lumps, each with 12 edges and 8 vertices. Two balls united at a touching point give 2 lumps of 1 face each, with no vertex at the touching point. A cylinder minus a ball whose equator touches the cylinder's wall gives 2 lumps, each with its own ring edge on that circle.
- One connected piece of material that touches itself only at isolated points is one lump; the points add no entity.
- One connected piece of material that touches itself along a curve, or at a point where its edges meet (the singular point of two meeting branches, such as the node of two cylinders tangent from inside), cannot be built as a manifold solid. It is refused with `boolean.self-touching-shell`, stage `sew`.

### Example: an exact answer

```json
{
  "format": "exact-brep-manifest/2",
  "question_id": "q-7f3a9c01",
  "kernel": { "name": "exact-brep", "version": "0.1.0", "commit": "0123abc" },
  "status": "answered",
  "exact": true,
  "result": { "kind": "solid", "lumps": 1, "shells": 1, "faces": 7, "edges": 14, "ring_edges": 2,
              "vertices": 8, "faces_by_kind": { "plane": 6, "cylinder": 1 },
              "edges_by_kind": { "circle": 2, "line": 12 } },
  "properties": {
    "volume": { "lo": "717.2566611769186085383", "hi": "717.2566611769186085384" },
    "area": { "lo": "...", "hi": "..." },
    "centroid": [ { "lo": "...", "hi": "..." }, { "lo": "...", "hi": "..." }, { "lo": "...", "hi": "..." } ]
  },
  "step": { "file": "answer.step", "max_deviation": "0", "approximated_curves": 0,
            "added_vertices": 2, "added_edges": 0 },
  "marks": [],
  "refusal": null,
  "stats": { "wall_ms": 3, "decisions": 412, "filter": 398, "refined": 0, "exact": 2,
             "coincident_exact": 12, "coincident_separation": 0, "coincident_fact": 0,
             "not_exact": 0, "refused": 0, "max_bits": 53, "work": 9120,
             "per_step": [ { "id": "s1", "wall_ms": 2, "decisions": 380, "not_exact": 0 } ] }
}
```

The ellipses stand for numbers the kernel computes; a real manifest never contains `"..."`.

### Example: a refusal

```json
{
  "format": "exact-brep-manifest/2",
  "question_id": "q-0c51e2d4",
  "kernel": { "name": "exact-brep", "version": "0.1.0", "commit": "0123abc" },
  "status": "refused",
  "exact": false,
  "result": null,
  "properties": null,
  "step": null,
  "marks": [],
  "refusal": { "code": "intersect.unsupported-pair", "step": "s1", "stage": "intersect",
               "detail": "cone and torus are not covered yet" },
  "stats": { "wall_ms": 0, "decisions": 6, "filter": 6, "refined": 0, "exact": 0,
             "coincident_exact": 0, "coincident_separation": 0, "coincident_fact": 0,
             "not_exact": 0, "refused": 0, "max_bits": 53, "work": 6,
             "per_step": [ { "id": "s1", "wall_ms": 0, "decisions": 6, "not_exact": 0 } ] }
}
```

## The STEP file

`answer.step` is written by piece 7: ISO 10303-21 text, schema AP242, millimetres, one product named after the question id. See `docs/pieces/07-step-export.md`.

- Each lump is one `MANIFOLD_SOLID_BREP`, or one `BREP_WITH_VOIDS` when it has voids, all in one `ADVANCED_BREP_SHAPE_REPRESENTATION`. Lumps that touch are written with their own faces, edges and vertices, so the file holds no non-manifold solid (ND4).
- Every face's surface is written as its analytic entity: `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE` or `TOROIDAL_SURFACE`. The file never holds a B-spline surface (ND7).
- Lines, circles and ellipses are written as `LINE`, `CIRCLE` and `ELLIPSE`, including those that arise where two surfaces meet. Only an intersection branch with no closed form is approximated, as one non-rational `B_SPLINE_CURVE_WITH_KNOTS` (ND5).
- Edges carry their 3D curve only. The file holds no `PCURVE`, `SURFACE_CURVE` or `SEAM_CURVE` (ND5). A pcurve would be a second approximation per face whose distance would need its own proof, and a reader can build it from the 3D curve and the analytic surface.
- `max_deviation` covers spline fitting and decimal rounding. It is at most `curve_tolerance` (default `1e-7` mm) unless the manifest holds a `step.tolerance-not-reached` mark, whose bound is then the proven deviation (ND5).
- The header holds no time of day, user name or machine name. Its `FILE_NAME` time stamp is always `1970-01-01T00:00:00`, so reruns give the same bytes.

## What the judge may rely on

- The files and exit codes above, and nothing else.
- Determinism, as stated.
- A refusal or a mark is an honest statement by the kernel: the kernel never returns a proven-looking answer (`exact: true`) that it did not prove.

## What the kernel may rely on

- The judge sends only questions in a format version listed by `exact-brep about`.
- The judge calls one `answer` per question in a fresh process.

## Decisions on the judge's open points

Each line gives the decision and where it is written. The dated entries in `INTERFACE-CHANGES.md` give the reasons.

| Point | Decision | Section |
| --- | --- | --- |
| ND1 | Rotations are quaternions `[w, x, y, z]` of exact numbers; integer quaternions are the usual form. Rotation matrices are not accepted. | Operands |
| ND2 | Questions use decimal (optional exponent) or `p/q` strings only. The manifest writes decimals only. | Units and numbers |
| ND3 | An empty answer has status `answered` with `result.kind` `empty`, properties 0, and no STEP file. It is never a refusal. | The manifest |
| ND4 | Touching pieces are separate lumps sharing no entity, each its own `MANIFOLD_SOLID_BREP`; a point contact inside one piece adds nothing; one piece touching itself along a curve, or at a point where its edges meet, is refused with `boolean.self-touching-shell`. | The answer's topology, The STEP file |
| ND5 | `max_deviation` is at most `curve_tolerance` (default `1e-7` mm) unless marked `step.tolerance-not-reached`; branches with no closed form are non-rational B-spline curves; no pcurves. | The STEP file |
| ND6 | Faces and edges are maximal; a merge that cannot be decided is refused. | The answer's topology |
| ND7 | Analytic surfaces and curves are always written as their analytic STEP entities. | The STEP file |
| ND8 | Questions carry no pose, family or cell label; `id` is opaque. | The question file |
| ND9 | Command line and exit codes as written: refusals exit 0, and every nonzero exit, timeout or missing manifest is a crash; 60 s and 4 GiB are the judge's outer limits; the kernel's default deadline is 50 s; empty sandbox; byte-identical reruns except `stats` and deadline refusals. | The program |
| ND10 | `properties` enclose the built body's exact properties, never the STEP copy's; they contain the true values unless a mark from outside pieces 6 and 7 is present. | The manifest |
| ND11 | Seams lie on the local half-plane `y = 0, x > 0`; sphere poles on the local `z` axis; torus tube seam on the outer equator. | Seams and poles |
