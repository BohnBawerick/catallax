# Numerical periodic cylinder admission

The retained analytic import path can join two complete cylindrical rims with
one explicit shared seam. The first source family has exactly two one-edge rim
loops, coaxial equal-radius circle carriers, opposite rim winding and existing
anchor vertices on one cylinder generatrix. Source vertices, rim edges and
analytic carriers remain unchanged. The derived line joins those anchors and
appears twice in one finite parameter loop, with opposite coedge senses and UV
branches separated by one period. A repeated cylindrical edge must pass this
periodic-cell check even if its UV width falls just below one turn; the ordinary
subperiod chart path cannot admit a failed shared seam. `NumericalSourceAdmission.topology_repairs`
records the source solid, face, rim edges and vertices, derived edge key and
versioned method. Entity numbers are source-local provenance, not replacement
correspondence.

A closed circle requires the same explicit endpoint vertex and the canonical
128-bit rational approximation of one turn. An open arc remains below one turn.
A full cylinder chart requires two closed rims and a single paired seam; an
arbitrary wide polygon remains refused. Existing finite winding, endpoint,
curve-on-surface samples and closed topology checks still run. A sloped paired
chart after decimal frame rounding is allowed only when both branches remain
translations by one period and its numerical residual checks pass. The STEP
writer validates the actual rounded chart before emitting paired uses as a
`SEAM_CURVE` with two cylindrical pcurves. Complete circular holes also require
a polygonal outer boundary and continuous disk-to-segment clearance above the
requested tolerance. The enclosing projected disk is conservative for tilted
carriers; safe cases can therefore refuse. A sampled rim polygon alone cannot
establish this clearance. Standalone complete circular caps remain supported;
other complete-circle multi-boundary arrangements refuse.

## Method and scope

This follows OCCT's [ShapeFix_Face::FixMissingSeam](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingAlgorithms/TKShHealing/ShapeFix/ShapeFix_Face.cxx)
recognition of opposite noncontractible wires and its preference for existing
edge endpoints when inserting a seam. The explicit one-edge/two-pcurve topology
matches [BRepLib_MakeFace](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingAlgorithms/TKTopAlgo/BRepLib/BRepLib_MakeFace.cxx).
The continuous boundary check follows OCCT
[BRepCheck_Face::IntersectWires](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Face.cxx)
in checking curve boundaries before classifying containment. This bounded
family uses a conservative analytic disk and interval-enclosed closest-segment
distance instead of adding a general curve intersector.
The narrower analytic family needs no general projection, new vertices or rim
splitting. P04's [SMLib consumer map](../P04/PORT-CONSUMERS.md) distinguishes the
unqualified `DropTracer` seam-start path from analytic cylinder operations; this
change does not enter that tracer or the spline Boolean/intersection paths.

Admission remains **numerical approximate**, at the requested positive tolerance
of at most 0.01 mm, with the unchanged 20M operation allowance. Sampled residuals
are not a complete body deviation bound or a global self-intersection proof.
All within-solid source face pairs retain their checked/unresolved coverage;
separate solids remain unexamined under the existing import policy. The inherited
optional curved proof uses a collective 1M prepayment; less remaining work can
still refuse rather than scale that reservation. Legacy planar predicate work
is not yet fully metered. This step does not claim general cylindrical healing,
spline healing, a surface-body result, a material Boolean or app replacement.

## Named reproduction and controls

Census source `00002345`, SHA-256
`8b138ab5f82850b19b017f86283b8074a88c7a271c4c4a01576a610bbfa7b40e`,
belongs to the original `self_intersecting` first-refusal cohort. The frozen
analytic census refused it as an unsupported numerical family. Independent
source inspection identifies two complete circular rim loops with compatible
anchors. Its geometry remains outside the repository.

Authored controls use square stock with a unit through bore. They check retained
source identity, one recorded derived seam, finite closed topology and STEP seam
output. Closed partial arcs, extra turns, implicit closure, independent seam
edges, wrong winding, double-period cells, near-turn shared seams, misaligned
anchors, between-sample circular crossings, near-feature clearance below the
requested tolerance, an over-cap source
gap and an exhausted request return no result. Existing analytic crossing and
source-pair isolation controls remain required. A full census after this repair
must keep the original seed, source hashes and denominator and distinguish
actual seam-repair admissions from other reader changes.

## Executed qualification

[qualification.json](qualification.json) records the current capture, its
**author-attested** association with kernel producer `51e9c0ffa13efcbcfb5c004b152e4c6014c6b5c1`, example executable hash,
independent checker hash and actual STEP bytes, including the authored bore hash.
The revision and executable-to-capture association come from the recorded
offload invocation and are author-attested for both the current and initial
named captures. The checker independently validates the supplied STEP artifacts,
compares response source hashes to fixed original pins, checks the named real
source file, and hashes the supplied executable. It does not independently
establish that this executable was built from the supplied revision or produced
those captures. Its own text hash identifies a checker version independently
of the attested kernel revision. This named-part provenance is distinct from
the census collector's actual Cargo-artifact build binding.
The earlier observed refusal is bound separately to
`dc0d7a4dab2c1e68feb52cb729800823ab4d97d3`. The named file imports as one numerical
body with 7 faces, 15 edges and 10 vertices. Import uses 1,190,064 work; source
STEP, pose rebuild and edited STEP use 219,540, 125,741 and 219,504 respectively,
totalling 1,754,849 under the unchanged 20M allowance. The
[initial named-part capture](initial-qualification.json) is retained with its own
producer, before the continuous-clearance and mandatory shared-seam guards.

The work-PC run passed 63 scoped release controls and 12 evidence controls. Both
actual exported files pass raw OCCT construction and whole-solid validity,
without healing, sewing or SameParameter. A source-carrier calculation accounts
for the source's metre units, checks the expected stock-minus-bore volume, and
verifies bore-void, material and outside probes before and after the pose. The
pose preserves volume and shifts the centroid by the requested translation. The
two earlier analytic real controls remain raw-valid with unchanged STEP hashes.
This qualifies a pose edit; it does not qualify arbitrary edits or a Catallax
Boolean on the part.

The authored bore's raw volume is `16-pi`, its bore is outside the material, and
its stock probe is inside. An OCCT reference cut produces a raw-valid result of
volume `8-pi/2` with the expected changed material probe. This is reference-kernel
Boolean evidence only. Run the retained checker against the actual private
captures with:

```sh
python tools/research/check_periodic_import.py \
  --source-revision 51e9c0ffa13efcbcfb5c004b152e4c6014c6b5c1 \
  --capture CAPTURE --imports IMPORTS --real-source 00002345.step \
  --example-binary numerical_import --out qualification.json
```

The [complete matching census](../../research/step-census.md) imports 400 of the
same 3,000 sources fully (13.33%), retaining all 378 prior full imports and adding
22 actual files with 44 recorded seams. That population count comes from all
3,000 terminal observations, separately from this named-part edit/STEP check.
The earlier prototype census remains private and contributes no current claim.

## Two seams in one real solid: STEP range repair

Census source `00003390`, SHA-256
`9ee5c0862a70c2f7eef35ce8bf61fb1a748fc173fa6c6b606745a9842e36e4a7`,
is one of the 400 admission-only results. It retains one solid with 16 faces,
40 source edges and 28 source vertices, plus two recorded derived seams.
The [before capture](multiseam-before.json) imports and exports under the frozen
`51e9c0ffa13efcbcfb5c004b152e4c6014c6b5c1` producer, but its source STEP is
raw-invalid: source face `#46` has a wire orientation result of
`BRepCheck_NotClosed`. A descending complete-circle parameter span was used for
its pcurve while the reader reconstructed an ascending natural full turn.
Before-fix pose STEP validity was not measured.

This repair follows OCCT's
[edge translation](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateEdge.cxx)
and [STEP edge writer](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/DataExchange/TKDESTEP/TopoDSToStep/TopoDSToStep_MakeStepEdge.cxx):
STEP edges do not carry their native parameter range, so vertex projection and
periodic range reconstruction must agree with their pcurves. The bounded
analytic writer puts a complete circle's reference phase at its existing shared
endpoint and uses the ascending natural turn. It retains the geometric center,
axis, radius and source vertex. Reversing a native edge direction also reverses
every coedge sense, preserving the physical boundary walk. The final serializer
uses that validated rounded definition. Partial circular arcs keep their
existing parameterization. This is an export repair; import admission and its
20M allowance are unchanged. It uses no general projection/tracer path from the
P04 SMLib consumer map above.

The [fixed capture](multiseam-qualification.json) pins kernel producer
`c9eb5861241fe67a7c3b8504ce62a588a9aa5a49`, executable/checker hashes and both actual
STEP hashes. Both exports pass raw OCCT validity without healing, sewing or
SameParameter. A separately labelled normal OCCT transfer of the actual source
may run healing; it supplies the independent volume/material reference. Source
and exported volume agree within relative `1e-8`; pose volume agrees within
relative `1e-9`, and centroid translation is checked within `1e-8` mm. Each
repaired bore has a void and stock probe before and after translation `(3,5,7)` mm.

The checker compares every source-local face, edge and vertex identity, all
source vertex positions within actual scalar rounding bounds, retained analytic
carriers, face-edge incidence, both derived seam endpoints and their opposite
uses. Complete-circle reference phases may rotate while their geometric
carriers and retained endpoints agree; partial circle reference frames remain
checked. A geometrically unchanged file with exchanged vertex identities is
rejected. These source-local checks do not establish correspondence between
changed source files. Public captures omit source coordinates and shape
properties; geometry stays in private proof storage.

Import, source STEP, pose rebuild and edited STEP use respectively 1,555,652,
604,772, 391,599 and 604,563 work, totalling 3,156,586 under the unchanged 20M
allowance. Linux formatting checks and 21 scoped native STEP tests passed,
including ten equivalent direction/period/phase variants. The nonzero-phase
authored bore also passes raw validity, volume and material checks. All three
earlier real workflows pass again in the separate
[post-fix regression capture](multiseam-regression-qualification.json); their
historical captures above remain unchanged.

The revision-to-executable association in these named captures remains
**author-attested**. The checker verifies the supplied executable hash, input
source hash, response source identity and actual STEP artifacts; it does not
independently establish build or run origin. No complete-body deviation bound,
arbitrary edit or Catallax Boolean is qualified. This named repair does not
supply an export-validity count for the 400 admissions. That count requires a
separate complete per-file measurement and will be published beside the
admission-only count.

```sh
python tools/research/check_multiseam_import.py \
  --source 00003390.step --capture PRIVATE_CAPTURE \
  --source-revision c9eb5861241fe67a7c3b8504ce62a588a9aa5a49 \
  --example-binary numerical_import \
  --example-sha256 ca103c7a8f352e1d4308b4cabd3027b823494e79a59a3bf1e4b70c2dcaa685cc \
  --phased-authored periodic-bore-phased.step --out qualification.json
```
