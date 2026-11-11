# R10 first step: what must be right for a proven body

Research date: 11 November 2026, against main
`54c4b1860ddca18397b4849cde2ea8629ad31cfb`. This is the audit that
[PLAN.md R10](PLAN.md#126-robustness-and-many-core-steps) places in M12's exit
run. It decides nothing by itself: the decision on the independent body checker
is taken after M12, from these numbers. No kernel behaviour, saved answer,
tolerance or work allowance changes. The only code is two read-only research
tools, [`tools/research/r10_census.py`](../tools/research/r10_census.py) and
[`crates/catallax-cli/examples/r10_audit.rs`](../crates/catallax-cli/examples/r10_audit.rs).

## 1. Answer

**Today the topology check re-derives none of the geometry a finished body's
validity rests on.** It confirms that each builder recorded the incidences it
built, and every face builder records them as it builds. On P0 to P5, 816 of the
872 incidence questions the check asks (94%) are answered by a fact some builder,
the Boolean or the intersection code recorded; the other 56 follow from a curve's
or point's own definition. With those facts removed, every one of the 816 fails.
The check never reads a face's sense and never compares two faces. Four planted
faults (an edge moved off its faces with its facts recorded, a face turned inside
out, two overlapping lumps, two coinciding spheres) all pass it.

The same facts also steer construction: 54,858 of the 126,155 exact-zero
decisions made while building P0 to P5 (43%) were settled by a recorded fact
rather than computed. A wrong fact is therefore not only passed by the check; it
can shape the body.

**So the code that must be right is every builder, not the numbers crate.**
For the three properties R10 names, that is about 111,000 code lines today:
`catallax-numbers`, `-shapes`, `-topology`, `-intersect`, `-boolean`, `-features`
and `-import`. The numbers crate is 14,204 of them (13%).

**A facts-blind checker is affordable for analytic bodies with code that exists.**
Rebuilt from the same geometry definitions with every recorded fact removed, the
14 recipes P0 to P5 attempt build 12 bodies (general cone and general torus
refuse). Every measured check passes on 7 of the 12. Edges on faces is established
on 11 (P5 refuses on its spline faces). Faces meeting only at edges is measured in
a narrower form, defined in section 4.2: the arrangement completes on 10 and
refuses on P3 and P5, and splits no face on any of the 10. The three bodies with
more than one shell are left undecided across shells: line tangent's two
cylinders meet along (5, 0, z), which the topology contract permits, and all
three have faces of different shells on one carrier, where a meeting over an
area would leave no section edge. Contacts inside one shell that split no face
are not covered. Inside and outside agrees on every face of 10 bodies; line
tangent has one disagreement at its contact (section 4.3) and P5 could not place
a point on five faces. The checks caught all four planted faults. Their cost ranges from 0.1
times the build's work (P1) to 25 times (crossing pipes, 99.3M units).

**Recommendation (section 6):** after M12, build the checker in three stages
that reuse the kernel's exact machinery without its facts: edges on faces first,
on every replay answer in CI; orientation second; face-against-face on request
and nightly. That takes the 40,301 lines of `catallax-features` and every
construction fact out of the trusted base. It does not reach "most of the numbers
crate": the reused numbers, shapes, intersection and Boolean code (61,499 lines)
stays common to the builder and the checker.

## 2. Method and reproduction

Two tools, both read-only.

`r10_census.py` counts non-test source: code lines per crate (not blank, not
only a comment), every site that forms a recorded fact (`FactOrigin::Construction`
or `FactOrigin::Proof`), every site that reads one, every site that sets a face
or coedge sense and every call of the topology check. Test files, benches,
examples and `#[cfg(test)]` modules are skipped. It is textual: a site is a line.

`r10_audit` builds each recipe through `catallax_features::evaluate`, as
`catallax-verify` does for a recipe, under the 100M replay allowance. For the
finished body it reports entity counts, the fact book by kind and origin, every
`Incident` question the topology check asks and which recorded fact answers it.
It then makes a **blind copy**: a fresh model with the same surfaces, curves and
points re-added in the same order, so every handle keeps its meaning, and the
result body's topology. The blind copy's fact book holds only what the store
records from a definition: an intersection curve on its two surfaces, a meeting
point on what it meets. The three properties are recomputed on the blind copy,
each pass under its own 100M allowance, and their work is reported in the
kernel's deterministic work units. A fault-injection section plants four faults
in 2 mm blocks and unit spheres. A `summary` block lists, by case name, which bodies completed and
which passed each check, so every denominator below can be read off the evidence.

```
python3 tools/research/r10_census.py > census.json
cargo build --release -p catallax-cli --example r10_audit
target/release/examples/r10_audit --workers 4 parts/P0/*.json \
  parts/P1-bracket.json parts/P1-drafted-bracket.json parts/P2-collar.json \
  parts/P3-filleted-boss.json parts/P4-shelled-tray.json parts/P5-flanged-housing.json
```

Cases run side by side, each with its own models and allowance. The report was
byte-identical at 1, 2 and 4 workers, SHA-256
`74bb0a207d15135040de2e6951b29553d5f91beb57a80e94b94a2bc1614eb143`; the run
takes 12 s on one worker. Both outputs are retained in [evidence](evidence/):
[`r10-2026-11-11-audit.json`](evidence/r10-2026-11-11-audit.json) and
[`r10-2026-11-11-census.json`](evidence/r10-2026-11-11-census.json). Raw run logs
are retained privately.

A work count of 0 means no charged decision was needed: the allowance charges
enclosure, sign and root work on exact reals, and those questions reduced to
rational identities. Wall-clock time is not a measure here.

## 3. (a) The code a proven body rests on

### 3.1 What the topology check reads and recomputes

`Model::check` ([`check.rs`](../crates/catallax-topology/src/check.rs), 310 code
lines, `check` at line 13) says of itself "Structural checks only. Incidence is
read from the fact book, never recomputed."

| The check | Reads | Recomputes |
| --- | --- | --- |
| Lumps, shells, faces, loops, coedges not shared; edges and vertices in one shell | Arena handles | Handle identity only |
| Loops closed | Coedge and edge senses, vertex handles | Handle equality at each join |
| Each edge used once in each sense per shell | Coedge senses | Counts |
| Shell edge-connected; Euler characteristic of the form 2 − 2g | Arena walk | Union-find, counts |
| Edge ends valid for the curve | `Curve::is_closed` | A flag of the curve kind |
| Edge curve on each face's surface; edge ends on the curve; vertex loops on the surface | `FactBook::holds(Incident, …)` (`check.rs:116`) | Nothing: a fact lookup with closure |
| Face sense, loop direction against the face normal, outward material | Not read | Not checked |
| Two faces meeting away from a shared edge | Not read | Not checked |

`FactBook::holds` ([`facts.rs:189`](../crates/catallax-numbers/src/facts.rs),
295 code lines) closes `Incident` over `Coincident` classes and follows chains of
incidence, so one wrong `Coincident` fact lends incidence to everything coincident
with its subject.

### 3.2 Who records the facts the check reads

Every face builder records the incidences of the face it builds as it builds
it, from the same data; only import proves them first:

| Recorder | Origin name | Where | What it asserts |
| --- | --- | --- | --- |
| Sweep face helper, used by `extrude.rs` and 31 other feature modules (revolve, draft, fillets, holes, chamfers, shells, lofts, polyhedra) | `features.sweep` | [`extrude.rs:49`](../crates/catallax-features/src/extrude.rs) `face`, recording at `:96` | Each edge curve on the new face's surface, each end on its curve, an apex on the surface |
| Primitive solids (block, cylinder, cone, sphere, torus) | `topology.<kind>` | [`build.rs:268`](../crates/catallax-topology/src/build.rs) `record`, `:278` `face` | The same, from exact rational sizes |
| Boolean sewing | `boolean.sew` | [`geometry.rs:164`](../crates/catallax-boolean/src/geometry.rs) `record_incident`, called from [`sewing.rs`](../crates/catallax-boolean/src/sewing.rs) lines 792, 813 and 829 | Every kept edge on every face that uses it, just before `checked` (`sewing.rs:998`) runs the topology check |
| Rigid copies | `topology.isometry` | [`transform.rs:32`](../crates/catallax-topology/src/transform.rs) | Copies each source fact to the image: a fact propagated, not recomputed |
| Surface pairs | `intersect.surfaces` | [`pair.rs:939`](../crates/catallax-intersect/src/pair.rs) | Closed-form conics and meeting points on both surfaces |
| Curve-surface hits | `intersect.curve-surface` | [`hits.rs:1289`](../crates/catallax-intersect/src/hits.rs) | Hit points on curve and surface |
| Store definitions | `shapes.definition` | [`store.rs:277`](../crates/catallax-shapes/src/store.rs) `add_curve`, `:296` `add_point` | An intersection curve on its two surfaces (true by definition); a meeting point on everything it meets, **including a `Position::Closed` coordinate that nothing checks against the objects it names** |
| Feature specials | `hole coaxial stations`, `draft.polyhedron-edge`, `fillet.cliff-foot`, `fillet.crossing-contact`, `features.curved-prism`, `topology.rational-atlas`, `boolean.bezier-interior`, `boolean.surface-probe` | `hole.rs:369`, `draft_build.rs:1535`, `fillet_build.rs:1694`, `fillet_equal_crossing.rs:364`, `end_condition/curved/build.rs:361`, `rational_atlas.rs:462`, `points.rs:991`, `probed_region.rs:46` | Family-specific incidence and coaxiality |
| STEP import | `import.step` | [`import/src/build.rs`](../crates/catallax-import/src/build.rs) lines 372 to 426 | Edge on face and vertex on edge, **proved exactly first** by `check_on_surface` (`:601`) and `check_on_curve` (`:673`) |
| Exact proofs | `proof: ExactArithmetic`, `proof: Algebraic` | `curve_pair.rs:32`, `graph.rs:124`, `imprinting.rs:67`, `points.rs:945`, `probed_region.rs:312`, `geometry.rs:183` | A relation decided at the time it was recorded |

The census finds 34 fact-forming sites in production code: 28 construction
sites and 6 proof sites. The STEP crate's 8 construction sites build export-only
chart models and do not touch a proven body. Facts are read at 46 sites, the
check's one among them; the rest are decisions that consult a fact before
computing: `GeomStore::side` (`store.rs:497`) and `point_on_curve` (`:541`),
`intersect_curve_surface`, which returns `Contained` when an `Incident` fact
holds (`hits.rs:1074`), the Boolean's region, probe and point-identity code.

One proof path launders a fact: `retain_incidence` (`curve_pair.rs:19`) and
`probe_at` (`points.rs:945`) record a `Proof` fact with the zero proof of the
decision behind it, and that decision may itself have been settled by a fact.
Crossing pipes' book holds 16 facts of origin `Proof(ConstructionFact)`.

Orientation is not in the fact book. It is the `sense` of each face, coedge and
edge, set at 28 face-construction sites in 20 files and 17 coedge sites in 13
files. The Boolean's `keep` stage (`stages.rs:268`) takes each kept piece's sense
from its operand and reverses a subtracted tool's.

### 3.3 Which recorded facts any independent check re-derives

| Independent check | Runs | Re-derives |
| --- | --- | --- |
| Sewing's shell orientation, `orient_shell` ([`sewing.rs:493`](../crates/catallax-boolean/src/sewing.rs)) | Every Boolean, unless its shells are already known to be outer | Sorts each new shell into outer or void by exact ray parity along one face's oriented normal. The senses of the shell's other faces are trusted through the check's edge-use rule. |
| Checked spline views: `SplinePrism::check`, `DraftedCylinder::check`, `RationalPatchFace`, `BezierExtrusionBand`, `LinearFilletBand`, `RollingMitreBand`, the rational atlas | When the verification path or an exporter meets a body of that family | Incidence, trims and material orientation for that family only (1,889 code lines) |
| Import's exact checks, including `check_solid`, a face-against-face self-intersection test for planar and some curved faces ([`import/src/intersect.rs:201`](../crates/catallax-import/src/intersect.rs), with `curved.rs`) | Each imported STEP solid | All three properties, from the file, before any fact is recorded |
| OpenCASCADE `BRepCheck_Analyzer` on the raw STEP file, no healing, 1e-7 precision ([`raw_step.py`](../tools/checks/raw_step.py)) | CI replay, for answers that export STEP | Curve-on-surface agreement within each edge's tolerance, wire and shell structure, in floating point, on the exported approximation rather than the exact body. As far as this audit read, it does not search for faces crossing each other; that would need OpenCASCADE's argument analyser, which is not run. |
| Field comparison | On request | Volume agreement with a kernel that never reads B-rep facts; a sampled whole-body check, not one of the three properties |

Everything else is re-derived by no check: the incidences recorded by
`features.sweep`, `topology.<kind>`, `boolean.sew`, `topology.isometry`,
`intersect.surfaces` and the feature specials; the closed-form meeting points of
`shapes.definition`; every face sense outside a Boolean's one face per shell; and
the absence of crossing faces in any body not imported from a file.

The recipe path of `catallax-verify` has no topology check of its own. The final
body is checked inside the step that made it when that step was a Boolean or one
of the 11 feature sites that call `Model::check`. The application path checks
every result (`catallax-app` `check_result`, `lib.rs:336`).

### 3.4 Measured on P0 to P5

| Case | Build work | Faces | Edges | Incidence questions | Answered by definitions alone | Zero decisions settled by a fact |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| block and sphere | 938,819 | 10 | 24 | 96 | 8 | 240 of 917 |
| crossing pipes | 3,897,430 | 20 | 20 | 56 | 8 | 399 of 2,117 |
| equal cylinder and sphere | 3,457 | 3 | 2 | 4 | 0 | 6 of 78 |
| line tangent | 131,052 | 6 | 4 | 8 | 0 | 68 of 308 |
| near miss | 4,202 | 6 | 4 | 8 | 0 | 2 of 78 |
| seam | 846,990 | 4 | 6 | 24 | 8 | 180 of 652 |
| P1 bracket | 42,152,303 | 19 | 42 | 144 | 0 | 25,088 of 56,777 |
| P1 drafted bracket | 42,377,532 | 20 | 44 | 152 | 0 | 25,354 of 57,039 |
| P2 collar | 2,751 | 6 | 6 | 12 | 0 | 0 of 54 |
| P3 filleted boss | 22,542 | 13 | 23 | 86 | 0 | 112 of 528 |
| P4 shelled tray | 1,554,975 | 15 | 36 | 144 | 32 | 2,819 of 5,423 |
| P5 flanged housing | 1,077,799 | 22 | 41 | 138 | 0 | 590 of 2,184 |
| **Total** | | | | **872** | **56** | **54,858 of 126,155** |

General cone and general torus refuse in the build (`intersect.traced-pair-required`)
and have no body to check. Every built body passes `Model::check`; its blind copy
fails it with exactly one `topology.missing-incidence` per question not answered
by a definition.

The direct origin of the 872 answers: `features.sweep` 382, `intersect.surfaces`
152, closure over other facts 102, `topology.block` 96, `shapes.definition` 56,
`topology.cylinder` 41, `boolean.sew` 40, `topology.isometry` 3.

### 3.5 Line counts

| Trusted for the three properties | Code lines today | With the staged checker (section 6) |
| --- | ---: | ---: |
| `catallax-numbers` | 14,204 | 14,204 |
| `catallax-shapes` | 26,914 | 26,914 |
| `catallax-intersect` | 8,897 | 8,897 |
| `catallax-boolean` | 11,484 | 11,484 (graph, imprint, region, containment and sewing reused) |
| `catallax-topology` | 3,904 | about 500 (model and entities); the check's incidence rule no longer matters |
| `catallax-features` | 40,301 | 0 |
| `catallax-import` | 5,569 | 0 |
| The checker | 0 | about 1,000, from the 781-line research example |
| **Total** | **111,273** | **about 63,000** |

`catallax-trace` (15,582 lines) returns unproved candidates only and no other
crate uses it yet; `catallax-measure` and `catallax-step` are trusted for
measurement and export answers, not for the body's validity.

## 4. (b) What a facts-blind checker recomputes, and its cost

Each property was recomputed on the blind copy. "Proved" means an exact outcome
agreeing with the body; "contradicted" means an exact outcome against it.

### 4.1 Every edge lies on both its faces

**Recompute:** for each coedge, the edge's curve lies on the face's surface; for
each edge, both ends lie on its curve; each vertex loop's point lies on its
surface. These are exactly the check's 872 questions, answered by computation.

**Reuse:** `catallax_intersect::intersect_curve_surface`, which answers
`Contained` when the curve lies in the surface; `GeomStore::point_on_curve` and
`GeomStore::side_of`, which compute when no fact holds. Store definitions answer
56 questions on their own.

**Not covered:** a curve proved to lie in the surface is stronger than the edge's
piece needing to; an intersection curve's branch is a further question the
equations do not settle. Neither occurred here.

### 4.2 No two faces meet except along shared edges

**The property:** for every pair of faces whose boxes overlap, where their
carriers meet, trimmed to both faces, must be their shared boundary and nothing
else. The topology contract does not ask this of faces in different shells:
"where material touches itself, each side keeps its own entities"
([pieces/03-topology.md](pieces/03-topology.md)), so two lumps may touch along a
curve with no shared edge. A checker must define which contacts are permitted.

**What was measured, in two parts, both narrower than the property:**

1. *No face split.* `catallax_boolean::arrange` with no candidate patches
   sections every face of the body against every other, imprints the meetings and
   hands each resulting cell to a classifier; keeping every cell re-sews the body.
   If two faces cross, at least one is cut into more than one cell. As many cells
   as faces shows that no meeting cuts a face; it does not show that every meeting
   lies on a shared boundary.
2. *No section and no shared carrier across shells.* Faces of different shells
   share no edge, so for them the property is that they do not meet at all. Each
   shell becomes a body of its own and the Boolean's `intersect` stage sections
   every pair of shells; any section edge or isolated point is a contact. An empty
   section graph is not enough on its own: where two faces lie on one carrier, the
   graph keeps only their boundaries, so two equal loopless spheres give an empty
   graph while coinciding everywhere. Every pair of faces from the two shells is
   therefore also tested for one carrier (the same surface, or
   `intersect_surfaces` answering `Coincident`), and any such pair leaves the
   shells undecided: whether the two faces' domains overlap is not measured. With
   no section and no shared carrier, faces of different shells do not meet, as far
   as the Boolean's section graph is complete for carriers that are not
   coincident. The part is vacuous for a body with one shell: 9 of the 12.

**Not covered:** a meeting inside one shell that splits no face, such as a face
touching another along a curve or at a point away from their shared edges.
`intersect_surfaces` alone is not enough for either part: it builds each meeting
curve's definition at no charged work, and the cost is all in trimming.

**Not covered:** carrier pairs the intersection engine cannot construct yet.
P3's arrangement refuses `intersect.traced-pair-required`. P5's refuses
`intersect.unsupported-pair`; of its 141 box-overlapping carrier pairs, 78 refuse
that way and 8 need tracing.

### 4.3 Inside and outside agree

**Recompute:** for every face, a point strictly inside it moved along the face's
outward normal (the surface gradient times the face sense) is outside the body,
and moved the other way is inside. This tests each face's sense against the
material, not only one face per shell.

**Reuse:** `catallax_boolean::face_point` and `classify_point`. The research
step is 10^-9 times the gradient. A checker would have to prove that no other
face lies within the step, which is property 2's work: on line tangent the
chosen point of one cylinder face is (5, 0, 5), on the line where the other
cylinder touches it, so the outward point is inside the other cylinder. That one
disagreement comes from the permitted contact the cross-shell section also finds,
not from a wrong face sense; it still means the property is not established there.

### 4.4 Cost in work units

| Case | Build | Edges on faces | Arrangement (cells / faces) | Across shells (sectioned; on one carrier) | Inside and outside | Check over build |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| block and sphere | 938,819 | 32,797 | 4,343,731 (10 / 10) | one shell | 576,191 | 5.3 |
| crossing pipes | 3,897,430 | 288 | 10,624,898 (20 / 20) | 2,967,394 (0; 8 face pairs) | 85,743,016 | 25.5 |
| equal cylinder and sphere | 3,457 | 0 | 2,444 (3 / 3) | one shell | 688 | 0.9 |
| line tangent | 131,052 | 0 | 131,620 (6 / 6) | 1,888 (1; 2 face pairs) | 1,848 (5 of 6, see 4.3) | 1.0 |
| near miss | 4,202 | 0 | 5,296 (6 / 6) | 1,022 (0; 2 face pairs) | 1,972 | 2.0 |
| seam | 846,990 | 19,331 | 1,221,566 (4 / 4) | one shell | 125,017 | 1.6 |
| P1 bracket | 42,152,303 | 1,080 | 4,055,329 (19 / 19) | one shell | 733,952 | 0.11 |
| P1 drafted bracket | 42,377,532 | 1,188 | 4,150,645 (20 / 20) | one shell | 737,146 | 0.12 |
| P2 collar | 2,751 | 0 | 8,488 (6 / 6) | one shell | 3,096 | 4.2 |
| P3 filleted boss | 22,542 | 720 | refused after 2,387: `intersect.traced-pair-required` | one shell | 844,247 | – |
| P4 shelled tray | 1,554,975 | 48,357 | 23,613,558 (15 / 15) | one shell | 7,248,255 | 19.9 |
| P5 flanged housing | 1,077,799 | 7,377,365, 20 refused | refused after 50,431: `intersect.unsupported-pair` | one shell | 24,363,991, 17 of 22 faces | – |

Each of the three bodies with two shells has one shell pair; the column gives
whether its section graph has an edge or point, and how many face pairs across
the two shells lie on one carrier.

Of the 12 built bodies, every measured check passes on 7: block and sphere,
equal cylinder and sphere, seam, both P1 brackets, P2 and P4, all single-shell
bodies. Crossing pipes and near miss are undecided across shells: their section
graphs are empty, but faces of the two shells share carriers (the void's
cylinders in crossing pipes, the cap planes of near miss's two cylinders), and
their domains were not compared. Line tangent has the permitted contact across
shells, two face pairs on one carrier and the orientation disagreement at the
contact. P3 refuses the arrangement. P5 refuses edges on faces on its
spline faces, refuses the arrangement and leaves five faces without an
orientation probe. Passing here means passing the measured checks, including the
narrower form of property 2; it does not establish that property in full.

The section across shells costs up to 3.0M units (crossing pipes, whose void is
a second shell). Edges on faces costs almost nothing on analytic bodies: at most 48,357 units,
under 3% of any build except P5. P5's 20 refusals are coedges on or along its
four spline faces (`intersect.parametric-surface` 16,
`intersect.parametric-curve` 4). Faces meeting only at edges is the expensive
property: up to 23.6M units, 15 times P4's build.
Inside and outside costs two point classifications per face, 85.7M on crossing
pipes; P5 could not place an interior point on five faces
(`shapes.bspline-domain` 4, `boolean.parametric-query` 1).

None of the three needed a fact a builder recorded. No allowance was exhausted.

### 4.5 The four planted faults

The first three were planted in a 2 mm block, with `Incident` facts recorded for
every question the check asks, as the face helpers do. The fourth is two unit
spheres from the primitive builder, which records its own facts.

| Fault | `Model::check` | Facts-blind recomputation |
| --- | --- | --- |
| The first edge's curve moved to a line 100 mm away | Passes | Edges on faces: 2 curves and 2 end points contradicted, 411 units |
| The first face's sense reversed | Passes | Inside and outside: that face contradicted, the other five refused `boolean.classify-ambiguous`; no face proved, 2,768 units |
| Two overlapping blocks as the two lumps of one body | Passes | Arrangement: 18 cells from 12 faces, then `boolean.self-touching-shell` at the sew, 94,983 units |
| Two equal spheres as the two lumps of one body | Passes | Across shells: empty section graph, one face pair on one carrier, so undecided rather than passed, 4 units |

## 5. Gaps in this audit

- The blind copy still trusts `shapes.definition`: a meeting point with a closed
  coordinate is incident to what it meets by declaration. A checker must also
  test those (`side_of` per surface), which is cheap.
- Property 2 is measured in a narrower form (section 4.2): a meeting inside one
  shell that splits no face is not detected; across shells, faces on one carrier
  leave the shells undecided because their domains are not compared, and contacts
  that are found are not judged against a definition of permitted touching.
- The orientation probe's step is not proved clear of other faces, and one P0
  case shows why it must be.
- The checker reuses the Boolean's arrangement and classification. A defect in
  that shared code could build a wrong body and also pass it. The field
  comparison, on request, is the only check with no shared code.
- P1 to P5 and P0 are small: 3 to 22 faces. Faces meeting only at edges grows
  with the number of box-overlapping face pairs (24 to 141 here), so larger parts
  will move its share.
- Only P0 to P5 were measured, not the 477 saved replay answers.

## 6. (c) Recommendation and trade-off

Build the checker after M12, as planned, in three stages, reusing the kernel's
exact machinery on a blind copy and never reading a recorded fact. Run it in CI
and on request, never on a click.

1. **Edges on faces, on every replay answer in CI.** Under 50,000 units on every
   analytic body measured; it caught the moved edge. It needs the closed-coordinate
   meeting points added (section 5) and a typed refusal kept for spline faces until
   their checked views can answer the same question without facts.
2. **Inside and outside.** Two classifications per face, at points proved clear
   of other faces. Up to 85.7M units on crossing pipes, so in CI it needs its own
   allowance; it caught the reversed face.
3. **Faces meet only at edges, on request and nightly.** First define which
   contacts are permitted: the topology contract lets lumps touch without sharing
   entities, as line tangent's do. The arrangement completes on 10 of the 12 built
   bodies, at up to 15 times the build's work, and caught the overlapping lumps,
   but shows only that no face is split. The checker still needs a step that
   matches every meeting curve and point against the shared boundary, which the
   arrangement's section graph contains but does not expose, and a comparison of
   face domains wherever faces share a carrier, without which every multi-shell
   body here stays undecided. Its refusals are pair families that M9, M10 and
   refusal repair already target.

The trade-off is the one the plan names, now measured. A full face-against-face
check is costly: it is the expensive property, its measured form refuses on 2
of the 12 built bodies today, and the full property needs a matching step not
measured here. Reusing the kernel keeps the check affordable and removes `catallax-features`,
the import builders and every construction fact (about 49,000 lines) from what
must be right for a proven body to be a valid solid. It leaves 61,499 lines of
numbers, shapes, intersection and Boolean code common to builder and checker.
Reaching "most of the numbers crate" would mean writing the checker's own face-pair
and classification code, as import did for planar and some curved faces in 630
lines; for general curved and spline faces that is a second Boolean, and not
recommended now.

Independent of R10, two findings could be acted on without the checker: the
recipe path of `catallax-verify` could run `Model::check` on every final body as
the application path does, and `Proof` facts could keep the origin of a
fact-settled decision rather than recording it as a proof. Neither is changed here.
