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

[qualification.json](qualification.json) retains the **initial named-part
capture, before the continuous-clearance and mandatory shared-seam guards**. It
binds that capture to kernel
producer `7e5b506f3d9448ff664fe4473c921eadc268bca5`, example executable hash,
independent checker hash and actual STEP bytes. The earlier observed refusal is
bound separately to `dc0d7a4dab2c1e68feb52cb729800823ab4d97d3`. The named file now
imports as one numerical body with 7 faces, 15 edges and 10 vertices. Import uses
1,189,708 work; source STEP, pose rebuild and edited STEP use 219,184, 125,385 and
219,148 respectively, totalling 1,753,425 under the unchanged 20M allowance.

The work-PC run passed 60 scoped release controls and 12 evidence controls. Both
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
  --source-revision 7e5b506f3d9448ff664fe4473c921eadc268bca5 \
  --capture CAPTURE --imports IMPORTS --real-source 00002345.step \
  --example-binary numerical_import --out qualification.json
```

The initial 3,000-file capture remains private because the stricter admission
guards require a fresh census. The complete matching rerun and later native
qualification belong in the subsequent census publication. No population
recovery count is inferred from this single qualified part.
