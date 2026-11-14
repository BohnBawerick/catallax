# FEA design for M16 and M20

**Status: agreed as written (12 November 2026).** The owner's 11 November
direction accepted the recommended D1-D4 choices for the first rectangular
cantilever below, and on 12 November 2026 the owner agreed this design as
written for the M16 groundwork and M20 solver lanes. D5 (benchmarks beyond this
cantilever) and D6 (imported parts) remain sequenced after the cantilever.
This document records the separate FEA design discussion required by
[PRD](PRD.md) and [PLAN](archive/PLAN-2026-11.md).

## First step: the field solves a rectangular cantilever

Build the simplest working field FEA demonstration first: a rectangular solid,
one fixed named end face and a defined load on the opposite named face. Solve
for its displacement, estimate **von Mises stress along the beam**, and compare
both stress and deflection with beam theory. Produce a stress picture and
numeric curves from the field solution. Beam-theory values are an independent
comparison, never a replacement for the field solve or its stress recovery.

Use this concrete starting configuration, retaining the previous draft's beam:

| Input | First cantilever |
| --- | --- |
| Geometry | Native rectangular beam `0 <= x <= 100 mm`, `-2.5 <= y,z <= 2.5 mm`; length `L = 100 mm`, width `b = 5 mm`, thickness `t = 5 mm`. |
| Material | One explicit homogeneous isotropic material, `E = 200000 N/mm²`, `nu = 0.3`. |
| Support | Fix all three displacement components on the surviving named face at `x = 0`. |
| Load | Uniform traction `(0, 0, -0.04) N/mm²` on the surviving named face at `x = L`, giving a total force `(0, 0, -1) N` over `25 mm²`. Other faces are traction-free; no body force or applied end couple. |
| Discretization | A Cartesian Q1 displacement grid aligned with the rectangle's faces, ordinary full-cell and face integration, and direct elimination of the fixed-face unknowns. Use the field definition and authoritative membership/cell queries to establish the material domain; do not substitute a B-rep volume mesh or a 1D beam element. |
| Output | Deflection along x and at the tip, the displacement-derived stress tensor, von Mises stress along x, load/support identities, reactions and actual numerical errors/work. Label the picture and numbers unproved numerical FEA. |

The accepted D1–D4 direction is sufficient for this aligned, single-material
case. Cut cells, sliver aggregation, curved-face Nitsche supports, general
material partitions, the complete six-query adapter, the plate and imported
parts are **not prerequisites for this first working cantilever**. They remain
later work below. Resolve the two named planar faces against the retained
field definition; arbitrary off-surface nearest-face and distance queries are
not needed to assemble this example. Do not guess when a required geometry
query refuses.

An aligned cell touching the surface can still receive `Unknown` from the
existing conservative closed-box classifier. For this native rectangle, a
separate exact check against its retained planar field constraints may certify
that the full cell is material up to its measure-zero boundary. Record and
charge that check; do not relabel the classifier's answer or infer full-cell
coverage from corner/centre samples. Without that certificate, retain the
unresolved outcome rather than dropping a boundary layer of the beam.

The separate [agreed cantilever implementation](M20-BEAM.md) runs the aligned
beam with displacement-derived tensors, incident-cell traces and a stress
picture. Its numerical errors and refusals remain visible; this is not a
completed M20 or a certified stress result. Refine the same beam and retain
each grid and its accuracy verdict. The implementation records
its grid, integration rule, solve policy, residual and resource limits before
running; it must not increase existing caller work allowances. Parallelism and
performance work are not prerequisites to the first result; whenever the
implementation uses 1, 2 and 4 workers, retain the required same-answer checks.

### Compare deflection and von Mises stress with beam theory

Recover strain and stress from the solved three-dimensional displacement,
`epsilon = sym(grad(u))`, `sigma = C(E, nu) : epsilon`. Compute von Mises from
that full stress tensor, including its shear components:

```text
sigma_vm = sqrt(((sigma_xx - sigma_yy)^2 + (sigma_yy - sigma_zz)^2
                + (sigma_zz - sigma_xx)^2) / 2
                + 3 * (tau_xy^2 + tau_yz^2 + tau_zx^2))
```

This is the equivalent-stress expression in the original
[CalculiX postprocessor documentation](https://web.mit.edu/calculix_v2.7/CalculiX/cgx_2.7/doc/cgx/node176.html).
It is numerical stress postprocessing, not a yield/safety qualification.

For the independent comparison use downward displacement magnitude and the
outer-fibre bending-stress magnitude:

```text
I                  = b*t^3/12 = 625/12 mm^4
w_EB(x)            = F*x^2*(3*L - x)/(6*E*I)
w_EB(L)            = F*L^3/(3*E*I) = 0.032 mm
abs(M(x))          = F*(L - x),  F = 1 N
abs(sigma_xx(x,z)) = F*(L - x)*abs(z)/I
sigma_vm,bend(x)   = F*(L - x)*(t/2)/I  at z = +/-t/2
```

The [FEAScript beam-method documentation](https://feascript.com/tutorials/cantilever-beam-1d.html)
gives the tip-loaded Euler–Bernoulli displacement curve. Roylance's original
[MIT stresses-in-beams notes, equations 5 and 7](https://web.mit.edu/course/3/3.11/www/modules/bstress.pdf)
derive rectangular-section inertia and the flexure stress. The last line above
is the uniaxial bending reference obtained by inserting that stress in the von
Mises expression; the field result must still use its full tensor. This beam
reference is an approximation to the loaded 3D solid, particularly near the
fixed and loaded faces.

Report area-mean downward tip displacement and the centreline deflection curve
`-u_z(x, 0, 0)`. For stress, publish the field-derived von Mises profile along
both outer-fibre paths `(x, 0, +/-2.5 mm)` together with the tensor components
and raw cell stresses. Evaluate the one-sided material trace of the displacement
gradient on each path; if multiple incident cells supply traces, report their
range instead of concealing it with colour smoothing. Show centreline von
Mises separately: bending normal stress vanishes there in beam theory, so it
cannot replace the outer-fibre comparison. The first bending comparison
stations are `x = 25, 50, 75 mm`, where the reference magnitudes are
`3.6, 2.4, 1.2 N/mm²`. They avoid the immediate end regions. Display end-region
stress as well, identifying it separately from these span comparisons; the
maximum at a clamped corner is not the beam-theory stress pass metric.

Keep signed and absolute deflection/stress errors at those declared locations,
reaction force/moment balance, residual and refinement history. A coarse first
cut can report the estimates and their actual errors; do not claim the existing
cantilever deflection mark until `abs(delta / 0.032 mm - 1) <= 0.02`. The
roadmap sets no separate cantilever von Mises percentage gate: report the
comparison and convergence without inventing a passed stress tolerance. Wrong
or unresolved results remain visible, with a typed refusal and counts when a
requested calculation cannot finish. General capability targets and the
477 saved replay answers are unchanged.

### What follows the cantilever

After this beam works, extend to cut-cell/sliver coverage and the wider M16
interface, then reopen D5's beyond-cantilever benchmarks and D6's imported-part
choice. M20 remains a pipeline proof, not a product feature or stress suitable
for signing off a part. Its later plate mark stays within 10% of the textbook
factor 3, and its imported-part run remains required for full M20 completion.
Deferral changes the implementation order, not those acceptance marks. The
stress picture must describe the part the B-rep exports. General wall-clock
gates wait for M26; field capability targets stay in their existing milestones.

## Source baseline and existing capability

The source inspection is of commit
`092ffab21a5fd168bbd881dc146326ae3cecbbc1`. The roadmap rewrite
([merged change](https://github.com/BohnBawerick/catallax/pull/361)) and the M14
hardware-preview decision
([merged change](https://github.com/BohnBawerick/catallax/pull/363)) are present.
The following are concrete APIs, rather than promises of M16 completion:

| Source | Existing contract and limit |
| --- | --- |
| [Field plan](../crates/catallax-field/src/lib.rs) | `Plan::proven_checked(&RPoint3, &mut Workspace, &mut Decider)` returns `Result<Query, Refusal>`. `Location` distinguishes `Inside`, `Outside`, `Boundary` and `Indeterminate`. `Value::Residual` is not Euclidean distance; `Query::distance` is a separate unsigned enclosure. `Direction` is a boundary normal only where certified, otherwise `Ambiguous` or `Undefined`. |
| [F4 interface](interface-changes/2026-11-11-unproved-processor-preview.md) and the field plan | `Plan::unproved_preview` returns `UnprovedValue`, labelled with `PreviewHardware::Processor`. It is display-only. `Plan::fast` and `fast_batch` still return tentative `Approximate` values; deterministic arithmetic alone does not make them proved. `true_distance_formula` is not an error certificate. No graphics-card evaluator is claimed here. |
| [Box classification](../crates/catallax-field/src/printing.rs) | Public `printing::classify_box(plan, lo, hi, dec)` classifies a closed rational box as `BoxState::{Inside, Outside, Unknown}`. `Unknown` means intervals did not settle it; it need not be a demonstrated geometric cut. The prepared `printing::Classifier` is crate-private. F6's reusable solver-facing preparation and target coverage must be delivered and documented, not inferred from this wrapper. |
| [Distance preparation](../crates/catallax-field/src/distance.rs) | `Plan::prepare_distances` returns `distance::Prepared`; `point_to_face`, `point_to_boundary` and `point_to_face_witness` work against surviving trims. Preparation needs a compact result bound and covered distance domains. A finite attained upper bound and convergence must be inspected; an unfinished search is not exact distance. |
| [F8 ownership](../crates/catallax-field/src/distance/nearest.rs), [explanation](FIELD.md#nearest-named-face) | `Prepared::nearest_face` returns `Ownership::Proved(Owner)` or `Ownership::Refused(Ambiguity)`, inside `Result`. `Owner` includes `origin`, witness, distance bounds, witness direction and `margin`. The margin is a lower bound on distance to an ownership border. It is not the exact border distance. The direction has no uniqueness or angular-error guarantee. |
| [Face registry](../crates/catallax-field/src/lib.rs) | `Face::origin` is a recipe name. `FaceId` is local to one immutable plan. `Plan::select` resolves a name; `face_query_checked` tests surviving trim membership and outward direction at a rational point. A registry entry alone does not establish that a boundary patch survives. |
| [Mass properties](../crates/catallax-field/src/mass.rs) | `mass::properties` computes certified whole-field volume, centroid and inertia with an estimate and work counts. It is not an exported per-grid-cell quadrature or material-region API. |
| [Field preview](../crates/catallax-field/src/display/preview.rs), [application preview](../crates/catallax-app/src/preview.rs) | `display::preview_plan` returns a distinct `PreviewMesh`. Application `preview_body`, `build_preview` and `construction_preview` expose display artifacts, not solver geometry certificates. |

F8 can refuse at equal-distance names, edges, uncovered carriers and unresolved
buried candidates. In particular, the documented buried-face limitation is not
fixed by this proposal. Current distance coverage includes planes, cylinders,
cones, spheres, ring tori and elliptic cylinders; spline carriers and horn or
spindle torus sheets can refuse with `shapes.distance-uncovered`. F4 sign checks
on finitely many points do not extend that coverage or certify solver inputs.
F7, material regions and the unified interface below are still requirements.
This read-only inspection supplies no new test, timing or benchmark results.

## Solver direction: beam first, cut-cell extension later

The scoped accepted D1 direction is static, small-strain, isotropic linear
elasticity on a uniform Cartesian background grid with continuous trilinear
(Q1) displacement elements and three components per unconstrained node. For
the first aligned rectangle there are no partially filled elements or cut-only
unknowns to aggregate. Here “no-mesh” means no body-fitted volume meshing;
the solver still has a grid. The later cut-cell extension integrates actual
material portions and uses aggregated basis constraints, rather than replacing
a partially filled element with a full one having scaled Young's modulus.

The proposed bulk problem is

```text
epsilon(u) = (grad(u) + transpose(grad(u))) / 2
sigma(u)   = C(E, nu) : epsilon(u)
sum_regions integral_region epsilon(v) : sigma(u) dV
           = integral_body v . b dV + integral_loaded_faces v . t dA
```

Prescribed displacements supply the support conditions. This small-strain weak
form and isotropic constitutive relation are described in the author's
[FEniCSx linear-elasticity implementation](https://bleyerj.github.io/comet-fenicsx/intro/linear_elasticity/linear_elasticity.html#variational-formulation).
The first beam needs one constant material, one face traction and one fixed
face only. Component constraints, material interfaces and other boundary types
can follow it; contact, plasticity, dynamics, buckling and nearly incompressible
materials require later designs. This is a design scope, not current capability.

For the later cut-cell extension, retain the D1 recommendation of aggregated
unfitted finite elements: constrain degrees of freedom supported only on cut
cells using an interior root cell, retaining
integration on every physical cut portion. Badia, Verdugo and Martín's
[original paper, sections 2–5](https://arxiv.org/html/1709.09122v1)
gives the aggregation and extension construction and conditional stability
analysis. Its numerical examples solve Poisson problems; its results do not
certify Catallax elasticity, curved integration, material interfaces or stress
recovery. Those need their own validation.

Use checked strong elimination for the first beam's grid-aligned fixed face.
Symmetric Nitsche enforcement of supports on cut named faces is deferred;
it is not a dependency of this first solve. De Prenter et al.'s
[finite-cell analysis](https://arxiv.org/abs/1601.05129)
documents small-volume conditioning, distinguishes small basis support from
near dependence, and studies preconditioning with Nitsche conditions in linear
elasticity. A diagonal rescaling alone is therefore not our assumed sliver fix.
Penalty selection and elasticity coercivity remain implementation evidence to
produce, not a theorem supplied by choosing the method's name.

Start with fixed-order assembly and the simplest inspectable linear solve that
fits the declared grid and existing allowances. Preconditioned conjugate
gradients and a separate direct-solve diagnostic are later options, not two
solvers that must exist before the first beam works. Record the chosen solve,
residual target, integration rule, iteration and cell limits before a run.
Loss of a valid system or exhaustion returns an explicit refusal and counts, never a plausible
successful picture. Q1 bending accuracy is a benchmark question, not a promise
that a coarse grid will pass.

## Later M16 interface: six queries

The following names describe contracts, **not existing Rust symbols**. The
accepted D2 recommendation is applied only to the first beam's needs now:
point/cell membership, full-cell and named-planar-face integration, stable face
resolution and a single explicit material assignment. Q4's general distance
route, cut moments/connectivity and multi-region Q6 are later work, not blockers
to the aligned cantilever. Full M16 should publish one prepared, immutable,
revision-bound adapter with per-worker scratch. Inputs include the retained
definition/import identity, exact rational grid origin and spacings, material
assignment revision, requested enclosure
widths and existing caller work policy. Preparation and queries charge that
policy without resetting or increasing it. Coordinates use millimetres; forces
use newtons and stresses/moduli use N/mm². Imported units are converted once
and retained in provenance.

Every answer identifies its revision, cell/point, units, evidence class, attained
bounds, convergence/coverage state and logical work. Preview suggestions live
in a different type. `Unknown`, boundary ambiguity, no admitted witness and
unsupported coverage remain distinct from success. A certified enclosure does
not assert that a solver's requested width was attained. Fixed query order and
canonical reductions preserve authoritative answers, refusals and work with
1, 2 and 4 workers; private caches must not make these history-dependent.

| Query | Required M16 contract | Existing starting point / missing work |
| --- | --- | --- |
| Q1: point membership | Exact point classification with `Inside`, `Outside`, `Boundary`, `Indeterminate`, plus typed numeric/coverage refusal. Boundary is never silently material or void. | `Plan::proven_checked` supplies the rational-point route. General algebraic boundary witnesses need an explicit adapter contract; do not round them into an exact membership claim. F4 can only supply a separately labelled display suggestion. |
| Q2: whole-cell state | For an exact closed axis-aligned box: certified all-material, certified all-void, or unresolved. Distinguish a proved cut with positive material and void measure from merely unresolved interval bounds. Shared grid facets have a canonical accounting rule; measure-zero contacts do not create volume or load twice. | `printing::classify_box` is concrete; prepared F6 reuse and coverage are groundwork. Corner or centre signs cannot prove the state of a whole cell. An unresolved box goes to Q3 or refuses within the fixed allowance. |
| Q3: cut-cell measure and integration | Enclose volume fraction `alpha` in `[0,1]`, volume and local first/second moments over material intersected with the cell, split by material region. Also provide the higher moments or integration data the selected element actually needs, and named-boundary integration data as specified below. Return achieved widths, unresolved measure and consumed work. | Whole-field `mass::properties` supplies useful integration precedent, not this per-cell API. F7 must deliver it. No clipping negative estimates, discarding unresolved pieces, or replacing integration with `alpha * full_cell_stiffness`. |
| Q4: distance and direction | Unsigned true-boundary distance enclosure, attained trimmed witness where available, convergence and direction semantics. Distinguish direction toward the witness, outward face normal and an actual certified unique nearest direction. Report ambiguity/undefined direction at ties, zero distance, edges and medial configurations. | `Query::distance`, `Prepared::point_to_boundary` and F8 witnesses provide pieces. `Query::direction` is a boundary normal, not an off-surface gradient. The fast true-distance-with-direction route required by M16 is not claimed to exist; residuals and F4 values cannot substitute. |
| Q5: owning named face | Nearest surviving recipe origin and trimmed witness, distance bounds, contenders on refusal, and a certified border-distance lower bound if available. Preserve ties; include name-to-result patch membership for boundary integration. | Reuse F8 `nearest_face`, including `Owner::margin` semantics and coverage refusals. Call `prepare_nearest` before sharing `Prepared` across workers. Exact distance to a border would be extra capability; M20 should accept a lower bound. |
| Q6: material region | At a point return void, one stable region/material id, interface membership, or unresolved/conflicting assignment. Over a cell return a region partition with Q3 integration for each portion and identified interfaces. Never infer material from nearest surface ownership. | No general region/material assignment API was found in the field plan. M16 must add one. Existing unit-density mass properties and lattice geometry are not a material model. |

### Integration data needed within Q3

Use cell-local normalized coordinates `xi = (x - cell_origin) / cell_spacing`.
Define `M_abc = integral_(cell intersect region) xi_x^a xi_y^b xi_z^c dV`.
`M_000` gives volume, with fraction `M_000 / cell_volume`; first and second
moments give centroids and inertia after the stated coordinate transformation.
They do not by themselves suffice for Q1 stiffness: for example, a product of
two x derivatives of Q1 shape functions includes `xi_y² xi_z²`.
This is a direct polynomial-degree check on the proposed element.

Recommend the full tensor set `0 <= a,b,c <= 2` (27 enclosed moments per region
per cell), or equivalently integration data with an explicit bound for those
polynomial integrals. Keep an estimate separate from each enclosure. Curved
cuts need bounded refinement, not an assumption that an ordinary Gauss rule on
the uncut cube integrates the cut domain. If widths cannot be met, report the
remaining volume/moment bounds and refuse the requested integration quality.
F7's volume acceptance remains: the sum of fraction times cell volume encloses
the B-rep volume on P1 and P4. A too-wide enclosure is not evidence of useful
solver accuracy; achieved widths must also be published. The planned 64³ grid
timing belongs to F7's own validation, not this document.

Q3 must identify disconnected material pieces and their certified facet
connections, or explicitly report that connectivity is unresolved. The later
cut-cell solver must refuse a cell/basis support that would couple separate physical
components across void. The proposed single Q1 space does not resolve such a
case just because its integrated volume is accurate. Refining a later request
may make that topology representable.

Q3 must also describe `cell intersect surviving named face`: covered patch,
area weights, outward normal data, boundary integration estimates and attained
error bounds sufficient for tractions and support terms. Q4/Q5 alone give
pointwise distance/ownership, not surface area integration. Surface normals can
vary on curved patches; normal-weighted integrals or a bounded quadrature rule
must support the chosen form. This is an explicit proposed addition to the
six-query contract, not a seventh query or a claim that F8 already integrates
faces. Integration on cell facets used to establish material connectivity must
likewise be available through Q3. The broader integration/connectivity scope
remains a proposal to agree before M16 can promise the full solver interface
is complete. D2's limited acceptance
for the first beam does not declare these capabilities implemented.

## After the aligned beam: cut cells and slivers (R14)

The updated M16 line asks the discussion to compare three approaches. Retain
cut-domain integration with aggregation as the D1 recommendation when extending
the beam to cut cells. D1 is accepted only as the first beam needs it; the
aligned case needs no sliver machinery. This section and its test do not block
that first result, and choosing a stabilising penalty instead would require a
later reviewed change.

| R14 approach | Assessment for M20 |
| --- | --- |
| Stiffness scaled by volume fraction | Reject as the default cut model. It retains only a zeroth moment of a varying stiffness integrand; Q3's polynomial example shows what it loses. A small fraction still carries real boundary force. |
| Finite cell method | Keep physical cut-domain integration as the starting principle, but require an explicit small-cut treatment; the cited de Prenter analysis documents conditioning concerns. Unaggregated finite cells with a dedicated preconditioner remain a D1 alternative. |
| Cut-cell method with a stabilising penalty | A viable alternative requiring a specified form and parameters. Burman's original [ghost-penalty paper](https://www.numdam.org/item/CRMATH_2010__348_21-22_1217_0/) describes penalty stabilization for fictitious-domain robustness. It does not establish that an arbitrary penalty or Catallax implementation passes the elasticity marks. |

Integrate interior cells normally, omit only certified void cells, and retain
every nonzero physical cut portion. `Unknown` is refined or reported, never
treated as void. Do not add fictitious material, clamp a small fraction upward,
erode a face to the grid, or discard a cell below a numerical fraction threshold.
Geometry error and the approximation in displacement are different quantities.

Recommend interior-root aggregation for **all** cut cells, avoiding a hidden
small-fraction cutoff. Constrain the basis; do not merge away the physical
volume, material or named boundary records. Use deterministic breadth-first
layers and a fixed global cell-id tie break to find roots. Connections must be
through a certified material-bearing shared facet of the same region, not
through a void gap, a different material, or a point/edge contact. Enforce
displacement continuity consistently at shared degrees of freedom, including
where different aggregates meet; independently replacing each cut cell's
displacements by a root's values is not the proposed finite element space.

Bound the permitted root path/aggregate extent in the request. If a thin
component has no eligible interior root, or connectivity cannot be proved,
refuse with the component, minimum fraction bounds, rootless-cell count and
work used. A later request can explicitly select a finer grid under unchanged
allowances. Do not fill the gap or silently raise resources. The root
existence/path assumptions in the aggregation paper matter; its theorem is
not a guarantee for every thin CAD feature. A ghost-penalty method is an owner
alternative (D1), not an automatic fallback that changes the accepted method.

### Sliver case required when cut-cell support is added

Use a manufactured affine-displacement elasticity patch on an exact rectangular
prism. With `h = 1 mm`, let its domain be
`[0, 4h + epsilon] x [0, 4h] x [0, 4h]` and the background x planes be integer
multiples of `h`. For `epsilon/h` in `{1/2, 1/16, 1/256, 1/4096}`, the layer
`x in [4h, 5h]` contains a fraction `epsilon/h` of material. In an interior y/z
cell its exact cut volume is `epsilon * h²`, and its outer x face has area
`h²`, however small the fraction. Repeat after a grid translation to prevent
success depending on one boundary alignment. All lengths exceed the modelling
tolerance; no geometric tolerance is changed to manufacture success.

For the manufactured reference choose dimensionless `a = 1/10000` and
`u = (a*x, -nu*a*y, -nu*a*z)`, giving constant uniaxial stress
`sigma_xx = E*a` and zero other stress components. One diagnostic prescribes
that affine displacement on all named faces. A second prescribes it only on
the named x=0 face, applies traction `(E*a, 0, 0)` on the named outer sliver
face, and leaves the four lateral faces traction-free. The analytic boundary
data belong to this test harness; they do not promise a general expression
language for product loads. Record the exact
volume/moments, boundary area and load resultant, root map, unconstrained and
constrained unknown counts, iteration count, residual, displacement/stress
errors, force balance and work for every fraction. The loaded sliver face
must not lose its force as `epsilon` shrinks.

Before running this later sliver test, declare numeric patch/residual/integration
thresholds in its reviewed test record. Those values are not set by the first
cantilever direction or by deferred D5. Passing requires the exact measures
inside their reported enclosures, nonzero volume retained, no lost or duplicated boundary resultant, the affine
solution within the agreed numerical target, and success or the same typed
refusal with 1, 2 and 4 workers. Report conditioning/iteration trends; do not
claim this finite sweep proves fraction-independent conditioning. As a
negative case, put a separate prism of thickness `h/4096` wholly inside a cell
with no interior root. It must return the proposed rootless-component refusal,
with numbers, rather than disappear. Proposed refusal categories here need
their own implementation/interface review; they are not existing error codes.

## Named loads, supports and materials

Store a boundary condition as recipe/import revision, stable origin selection,
condition kind, components, coordinate frame and values with units. Resolve
all surviving patches carrying the selected name at that revision; do not save
a registry index or a preview triangle as the persistent selection. A split
face retains its named pieces; an absent, ambiguous or unsupported selection
must be reported, not rebound to the nearest convenient face after an edit.
Preview picking may suggest a name, but authoritative selection validation
starts independently from the retained definitions.

Integrate traction over Q3's selected boundary patches, including cut faces.
A prescribed total force needs an explicit distribution rule; recommend
uniform traction normalized by the selected face's area, with the normalization
error recorded. Pressure would require the outward face normal, not F8's
query-to-witness vector. Report assembled force and moment about a recorded
origin and support reactions. Do not smear a load through a band of voxels or
apply a condition to every nearby differently named face. Boundary edges have
zero area in these face loads; ownership ties there do not justify dropping a
finite unresolved patch. Incompatible component constraints on overlapping
selections refuse; the treatment of compatible duplicates is explicit and
idempotent. Validate removal of rigid modes before solving; a singular system
is a refusal, not zero stress.

D4 is accepted for the first beam's explicit whole-part homogeneous material
and grid-aligned fixed support only: `E = 200000 N/mm²`, `nu = 0.3`. No region
partition or two-material continuity test is a prerequisite for this beam.
Full M16 still needs region identities, assignments and interface queries;
using a single material now does not claim that general Q6 is complete.

For later multi-material work retain the proposed disjoint, piecewise-constant,
perfectly bonded regions, with `E > 0` and `-1 < nu < 1/2`; sum their separate
stiffness integrals with shared displacement. This does not promise accurate
nearly incompressible Q1 behaviour. Region gaps, overlaps and unknown labels
refuse; do not average moduli across an interface. Aggregation must respect
region boundaries. A two-region integration/continuity check is required before
claiming that later capability. Region/Boolean assignment semantics and curved
support enforcement remain outside D4's limited acceptance. Contact and
debonding require another design.

## Deferred benchmarks and imported parts: D5 and D6

**Reopen these only after the native cantilever works.** The cantilever's
configuration, stress/deflection comparison and unchanged 2% deflection mark
are specified first above; D5 does not postpone that work or make a plate run
its prerequisite. The following retained configurations are proposals for
later discussion, not accepted first-step requirements or measurements.

| Deferred case | Retained proposal | Unchanged full-M20 requirement |
| --- | --- | --- |
| D5: plate with a hole | Thin 3D plate `200 x 200 x 1 mm`, radius `5 mm`; homogeneous material as above, free thickness faces. Quarter model with named x=0/y=0 component supports, remote x traction `sigma_inf = 1 N/mm²`, free remote y face and an explicit out-of-plane datum. Publish the hole's mid-thickness circumferential-stress profile and raw cell stresses. Freeze the stress evaluation rule, grid/refinement series and finite-domain comparison after the cantilever works. | `Kt = sigma_theta_theta,max / sigma_inf`; `abs(Kt / 3 - 1) <= 0.10`, i.e. `2.7 <= Kt <= 3.3`. A net-section nominal stress or colour maximum does not replace this definition. |
| D6: imported part | Repeat an accepted native configuration on the same imported geometry, with the same material, grid and named conditions. Retain analytic STEP via M21 as the recommendation, closed planar mesh via M22 as an alternative, and both as the broader option. Choose after the native cantilever works. | One complete imported solver run, not merely successful import, plus the required scoped cross-kernel agreement on the same part. |

The official [COMSOL Kirsch benchmark](https://doc.comsol.com/6.3/doc/com.comsol.help.models.sme.kirsch_plate/kirsch_plate.html)
uses a plane-stress model and distinguishes a finite plate from an infinite
domain. Its vertical-symmetry-line equation gives
`sigma_x(a) = sigma_inf/2 * (2 + 1 + 3) = 3*sigma_inf` at hole radius `a`.
The later 3D plate must investigate thin-plate and finite-domain error; the
reference does not establish that those errors vanish.

For those later cases, freeze dimensions, boundary selections, integration and
solve limits before results exist. Publish grid refinements and shifted-grid
cuts, force/moment balance and actual errors. Compare a larger plate domain at
matched near-hole resolution. No tolerance, reference or work allowance is
adjusted after seeing a picture. If a required query refuses, retain the typed
refusal, numbers and unfinished milestone requirement.

The imported run must retain source units, face names and import certificates.
A successful import alone does not establish M16 query coverage. At that later
stage bind the B-rep qualification/comparison and numerical solver report to
the identical definition/import revision; a picture of an older displayed
revision does not establish agreement. Do not promote a sampled comparison to
whole-solid proof. This deferral neither changes import tolerances nor claims
full M20 complete when only the native cantilever has run.

## Proof boundary and reporting

Follow the existing [M14 decision](ARCHITECTURE.md#unproved-hardware-preview-evaluator-decision-of-2026-11-10).
F4 hardware values and preview caches cannot choose proved membership, topology,
trim, ownership, integration acceptance or proof work ordering. They are not
certified FEA inputs. The scoped accepted D3 choice keeps them entirely out of
the first beam's assembly; retain that boundary for later M20 work.
They may colour an explicitly unproved display whose status is retained.

| Evidence | Meaning |
| --- | --- |
| Proven field geometry | Only the claims actually established by deterministic query enclosures, exact definitions and scoped certificates. A refused or unresolved query stays so. |
| Numerical solver output | Displacements and stresses from a declared discretization, integration approximation and algebraic solve. Even with certified geometry, they are unproved numerical results. A residual measures the discrete solve; it is not an error bound against continuum elasticity. |
| Benchmark agreement | Finite validation of the declared pipeline and configurations. The 2%/10% marks are not a theorem or a safety qualification for arbitrary parts. |
| Stress picture | Label `unproved — M20 pipeline demonstration`, with revision, material/condition identity, grid and current refusal/convergence status. No export qualification or proven-measurement slot receives solver values. |

Keep every proven answer and the 477 saved replay answers unchanged. Do not
relax modelling/approximation/import tolerances or the existing work allowances.
Geometry preparation and queries must reproduce answers, refusals, evidence
ordering and logical work with 1, 2 and 4 workers. Recommend a canonical solver
assembly/reduction order as well; do not claim that M14's display exception
already authorizes nondeterministic numerical FEA. Any broader arithmetic
exception needs its own owner ruling (D3). No clock chooses an answer.

The future run record should include immutable input hashes, geometry proof
state, imported source identity where applicable, grid and region partition,
boundary origins, force/moment checks, integration widths, aggregation choices,
linear-solver policy, errors and work. A refusal includes its typed reason,
attained bounds, counts and remaining requirement. Proposed reasons include
rootless component, incomplete boundary integration, invalid material partition
and unconstrained rigid mode; no error code is added by this proposal.

## Scoped decisions and deferred choices

The owner accepted the recommended D1–D4 options **only as needed for the first
rectangular cantilever**. Record that supplied direction without extending it
to the whole solver design or closing shared decision records. D5 and D6 are
deferred until the cantilever works. The designated reviewer must review this
revision before merge; neither CI nor this document substitutes for that review.

| Id | First-step status | Later choices and recommendation |
| --- | --- | --- |
| D1 | Scoped acceptance: Q1 Cartesian linear elasticity. The first grid-aligned rectangle has full cells, so no aggregation or stabilising penalty is needed for it. | After the beam works: (A) cut integration with interior-root aggregation; (B) ghost-penalty cut elements; (C) unaggregated finite cells with a dedicated small-cut preconditioner. Retain A as the recommendation with R14 sliver validation; no deletion of tiny cells or fictitious stiffness. Limited acceptance does not validate a general cut-cell implementation. |
| D2 | Scoped acceptance: use the proposed integration contract as far as the beam needs membership, full-cell/planar-face integration and named boundary resolution. Do not wait for all six general queries. | Later: (A) extend Q3 to higher moments, named curved-boundary integration and material-piece/facet connectivity; (B) a separate reviewed integration service; (C) remain aligned-only. Recommend A for full M16; C is sufficient only for the first beam and cannot discharge the later curved plate. Broader interface details remain proposals. |
| D3 | Scoped acceptance: deterministic authoritative geometry inputs and canonically ordered numerical assembly/solve, visibly unproved. Keep F4 preview arithmetic out of the beam's assembly. | Later alternatives are hardware-varying numerical assembly/solve or separately labelled approximate-geometry exploration. Recommend retaining the current boundary; either broader exception needs its own written ruling. No continuum-stress certification is granted. |
| D4 | Scoped acceptance: one explicit homogeneous isotropic material, one fixed named face with checked aligned elimination and one defined traction on a named face. | Later: (A) disjoint perfectly bonded piecewise-constant regions plus cut-face supports; (B) continue homogeneous-only while Q6 groundwork proceeds; (C) enlarge to contact/discontinuous displacement or another support formulation. Retain A as the broader recommendation; region/Boolean semantics, curved support enforcement and multi-region validation are not prerequisites for the beam and are not approved by this limited direction. |
| D5 | **Deferred until the cantilever works:** benchmarks beyond that cantilever. Its own deflection and von Mises comparison proceeds now. | Then discuss (A) the retained 3D plate configuration; (B) a 2D plane-stress Kirsch run plus a separate 3D demonstration; (C) other recorded configurations preserving the 10% factor-3 mark. Retain A as the recommendation. Freeze its stress rule, numerical targets and finite-domain check at that time; do not invent passed thresholds. |
| D6 | **Deferred until the native cantilever works:** imported-part selection and solver demonstration. | Then choose (A) analytic STEP via M21; (B) closed planar mesh via M22; (C) both. Retain A as the recommendation, B as the alternative with full required query coverage. One complete imported run remains the full-M20 minimum; nothing imported is required before the first native beam result. |

Implement and assess the field cantilever first. After it works, extend the
agreed scope through M16's cut-cell, named-boundary and material groundwork,
then the deferred native plate and imported evidence. Product FEA, trusted
stress error bounds and performance optimisation remain later work.
