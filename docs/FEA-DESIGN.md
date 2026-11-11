# FEA design for M16 and M20

**Status: proposal awaiting owner agreement.** This is the draft for the
separate FEA design discussion required by [PRD](PRD.md) and [PLAN](PLAN.md).
Publishing or merging this document does not record agreement with its choices.
Explicit design-discussion approval is required before this proposal is merged;
M16 cannot finish until the design is agreed. The decisions at the end remain
open. No solver or new field interface is implemented by this document.

M20 proves the pipeline from recipe to stress picture. It is not a product
feature or a stress result suitable for signing off a part. Its fixed marks are
cantilever deflection within 2% of the beam formula, plate-with-hole peak stress
within 10% of the factor 3, loads and supports on named faces, and the same run
on an imported part. The picture must describe the part the B-rep exports.
Capability comes first. General wall-clock gates wait for M26; the field
capability targets already assigned to F4, F6 and F7 stay in those milestones.

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

## Proposed solver family

Recommend static, small-strain, isotropic linear elasticity on a uniform
Cartesian background grid: continuous trilinear displacement elements (Q1
hexahedra), integrated over the actual material portion of each cell. There are
three displacement components per unconstrained node. Use cut-cell integration
and aggregated basis constraints at the boundary, rather than treating a
partially filled voxel as a full cell with a scaled Young's modulus. A body-fitted
tetrahedral mesh is not an M20 dependency. Here “no-mesh” means no body-fitted
volume meshing; the solver still has a grid and integration partitions.

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
M20 would initially support constant region properties, face tractions, fixed
supports and component constraints; contact, plasticity, dynamics, buckling and
nearly incompressible materials require later designs. These are scope
recommendations, not current capabilities.

For slivers, recommend aggregated unfitted finite elements: constrain degrees
of freedom supported only on cut cells using an interior root cell, retaining
integration on every physical cut portion. Badia, Verdugo and Martín's
[original paper, sections 2–5](https://arxiv.org/html/1709.09122v1)
gives the aggregation and extension construction and conditional stability
analysis. Its numerical examples solve Poisson problems; its results do not
certify Catallax elasticity, curved integration, material interfaces or stress
recovery. Those need their own validation.

Propose symmetric Nitsche enforcement of supports on cut named faces, with
strong elimination allowed for exactly grid-aligned supports as a checked
special case. De Prenter et al.'s
[finite-cell analysis](https://arxiv.org/abs/1601.05129)
documents small-volume conditioning, distinguishes small basis support from
near dependence, and studies preconditioning with Nitsche conditions in linear
elasticity. A diagonal rescaling alone is therefore not our assumed sliver fix.
Penalty selection and elasticity coercivity remain implementation evidence to
produce, not a theorem supplied by choosing the method's name.

Start with a fixed-order sparse assembly and an inspectable linear solve.
Propose preconditioned conjugate gradients after symmetry and positive
definiteness checks, with a small direct-solve reference for diagnosis. Solver
choice, residual target, integration targets, iteration and cell limits must be
recorded before benchmark runs; they are not selected here. Loss of a valid
system or exhaustion returns an explicit refusal and counts, never a plausible
successful picture. Q1 bending accuracy is a benchmark question, not a promise
that a coarse grid will pass.

## Proposed M16 interface: six queries

The following names describe contracts, **not existing Rust symbols**. M16
should publish one prepared, immutable, revision-bound adapter with per-worker
scratch. Inputs include the retained definition/import identity, exact rational
grid origin and spacings, material assignment revision, requested enclosure
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
connections, or explicitly report that connectivity is unresolved. The first
solver must refuse a cell/basis support that would couple separate physical
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
likewise be available through Q3. The owner must agree this scope before M16
can promise the solver interface is complete (decision D2).

## Cut cells and slivers: R14 proposal

The updated M16 line asks the discussion to compare three approaches. This
proposal recommends cut-domain integration with aggregation; the extra choice
of algebraic aggregation versus a stabilising penalty remains D1.

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

### Required sliver case

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

Before implementation runs, agree numeric patch/residual/integration thresholds
(D5). Passing requires the exact measures inside their reported enclosures,
nonzero volume retained, no lost or duplicated boundary resultant, the affine
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

Recommend an explicit whole-part material assignment for the two M20 marks,
with `E > 0` and `-1 < nu < 1/2`; benchmark `nu = 0.3`. This does not promise
accurate nearly incompressible Q1 behaviour. M16 still needs region identities,
assignments and interface queries: whole-part default material is an explicit
assignment, not an excuse to omit Q6. Recommend disjoint, piecewise-constant,
perfectly bonded regions; sum their separate stiffness integrals with shared
displacement. Region gaps, overlaps and unknown labels refuse. Never average
two moduli because a cell straddles an interface. Aggregation must respect
region boundaries. Require a two-region integration and continuity check
before claiming material-region support; failure to find a same-region root
remains explicit. Contact, debonding and discontinuous displacement are outside
this proposed M20 scope. The region schema and Boolean assignment semantics
are for owner decision D4.

## Benchmarks and acceptance

These are proposed reproducible configurations for the owner's agreement,
not measurements. Dimensions, grid origins/refinements, material, boundary
selections, integration targets and solve limits must be frozen before running
them. Keep both signed errors and absolute errors, all refusals and actual
counts. A failed mark stays failed; no reference, tolerance or allowance is
adjusted after seeing the stress picture.

| Case | Proposed setup and observable | Fixed roadmap pass mark |
| --- | --- | --- |
| Cantilever | A `100 x 5 x 5 mm` rectangular beam, length along x; `E = 200000 N/mm²`, `nu = 0.3`. Fix all displacement on named x=0 face. Apply a uniform transverse traction on named x=L face with resultant `F = 1 N`, no artificial end moment. Measure area-mean transverse tip displacement on that face, and report centreline displacement separately. Use `I = b*t³/12`, so `delta_EB = F*L³/(3*E*I) = 0.032 mm`. | `abs(delta / delta_EB - 1) <= 0.02`. Small-strain 3D end/shear effects and grid error count against this mark. Clamp stresses are not the benchmark observable. |
| Plate with a hole | Proposed thin 3D plate, full dimensions `200 x 200 x 1 mm`, hole radius `5 mm`, homogeneous material as above, free thickness faces. Use a quarter model with named x=0 and y=0 component supports, remote x traction `sigma_inf = 1 N/mm²`, free remote y face, and one explicit out-of-plane datum to remove translation. Resolve the circular hole with cut integration. Measure the maximum circumferential tensile stress on the hole boundary at mid-thickness using a fixed one-sided material stress evaluation rule; publish the angular profile and raw cell stresses. | `Kt = sigma_theta_theta,max / sigma_inf`; `abs(Kt / 3 - 1) <= 0.10`, i.e. `2.7 <= Kt <= 3.3`. Neither net-section nominal stress nor a colour-scale maximum replaces this definition. |

The beam reference is the tip-load Euler–Bernoulli solution documented by the
[FEAScript method tutorial](https://feascript.com/tutorials/cantilever-beam-1d.html).
Using a face resultant and a slender 3D solid is our proposed adaptation, not
an assertion that beam theory is exact at a finite 3D clamp.
The official [COMSOL Kirsch benchmark](https://doc.comsol.com/6.3/doc/com.comsol.help.models.sme.kirsch_plate/kirsch_plate.html)
uses a plane-stress model and distinguishes a finite plate from an infinite
domain. Its equation for stress on the vertical symmetry line gives
`sigma_x(a) = sigma_inf/2 * (2 + 1 + 3) = 3*sigma_inf` at the hole radius `a`.
The proposed 3D plate must demonstrate its thin-plate and finite-domain error;
this reference is not proof that those errors vanish.

For each case run at least three uniform grid refinements, including shifted
grid origins that create genuine cuts, and report unknowns, fractions, aggregate
sizes, integration widths, equilibrium residual/reactions and the observable.
Use the finest declared result for the mark, alongside its convergence history;
do not choose whichever coarse grid happened to cross the reference. For the
plate also compare a larger remote boundary at matched near-hole resolution.
Agree the stress evaluation rule, domain-size comparison target and refinement
series in D5 before results exist. Colour smoothing never supplies the peak.
Failure of Q1 to meet either fixed mark is reported and discussed, rather than
quietly introducing a beam/shell solver or changing the target.

Repeat at least one complete accepted benchmark on an imported representation
of that same geometry, with the same material, grid and named conditions. Prefer
the analytic STEP beam under M21; a validated closed planar mesh of the beam
under M22 is the alternative (D6). Build both kernels from the imported source,
retain source units, face names and import certificates, and request the
existing scoped cross-kernel agreement. Compare geometry, selections, load
resultants and the benchmark observable as well as the picture. A successful
import alone does not establish M16 distance/ownership/integration coverage;
if any required query refuses, report it and leave the imported M20 gate open.

The B-rep qualification/comparison and solver report bind to the identical
definition/import revision. A picture of an older displayed revision cannot
stand in for that agreement. No sampled comparison is promoted to whole-solid
proof. The import cap, any recorded override and local certificates remain
visible; this design creates no new import tolerance.

## Proof boundary and reporting

Follow the existing [M14 decision](ARCHITECTURE.md#unproved-hardware-preview-evaluator-decision-of-2026-11-10).
F4 hardware values and preview caches cannot choose proved membership, topology,
trim, ownership, integration acceptance or proof work ordering. They are not
certified FEA inputs. Recommend keeping them entirely out of M20 assembly (D3).
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

## Open owner decisions

All entries are unresolved. Recommendations support the discussion; they do
not record owner approval or close shared decisions.

| Id | Question and options | Recommendation and consequence |
| --- | --- | --- |
| D1 | Solver/slivers: (A) Q1 Cartesian elasticity with interior-root aggregation; (B) cut finite elements with ghost-penalty stabilization; (C) unaggregated finite cells with a dedicated small-cut preconditioner. | A, retaining every physical cut portion, plus the R14 test and explicit rootless refusal. B/C require a separate stabilization/preconditioner design and evidence before replacing A. No fraction deletion or fictitious stiffness fallback. |
| D2 | Are six point/cell queries sufficient as currently phrased? (A) extend Q3 to bounded higher moments, named-boundary integration and material-piece/facet connectivity; (B) add a separate reviewed integration interface; (C) restrict the first solver to grid-aligned boundaries. | A keeps one six-query adapter while stating what assembly actually needs. B is viable if the owner prefers a separate service. C cannot discharge the plate's curved-cut pipeline requirement alone. Agree this before M16 completes. |
| D3 | Arithmetic/proof policy: (A) certified deterministic geometry inputs and canonically ordered numerical solver output, visibly unproved; (B) allow hardware floating-point variation in numerical assembly/solve under a new explicit exception; (C) also permit approximate geometry in exploratory FEA under a separately labelled contract. | A for M20. Keep F4 out of assembly and preserve all existing proven answers. B/C need a new written boundary; M14's display permission does not decide them. No claim of continuum stress certification under any option. |
| D4 | Materials and supports: (A) piecewise-constant disjoint bonded regions with explicit whole-part assignment, named-face Nitsche supports and checked aligned elimination; (B) homogeneous-only solver initially, but finish Q6 groundwork; (C) enlarge scope to contact/discontinuous materials or another support formulation. | A as the target contract; validate each region separately and refuse absent roots/conflicts. B can stage M20's homogeneous marks without calling material-region groundwork complete. C needs another design. Owner must agree region/Boolean assignment semantics and support formulation. |
| D5 | Benchmark freeze: (A) the 3D beam/plate configurations above, with a declared stress rule and convergence/integration/solve targets; (B) a 2D plane-stress Kirsch run plus a separate 3D imported demonstration; (C) other documented sizes/BCs that preserve the fixed 2% and 10% marks. | A exercises the same six-query 3D path throughout. Freeze grid/refinement series, numerical thresholds, Nitsche policy and finite-domain comparison target in the discussion before implementation runs. Those numerical choices remain unset here; no benchmark has passed. |
| D6 | Imported gate: (A) repeat the beam as analytic STEP via M21; (B) repeat it as a closed planar mesh via M22; (C) require both. | A, with B acceptable if its full M16 query coverage lands first; record the selected path before acceptance. One complete imported run is the roadmap minimum; import success without solver queries is insufficient. |

After agreement, M16 should publish the reviewed interface and its coverage and
refusal tests, including F7's volume enclosure, named-boundary integration,
material partition and R14 slivers. M20 then implements the selected solver and
records the native and imported benchmark evidence. Product FEA, trusted stress
error bounds and performance optimisation remain later work.
