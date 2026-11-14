# M25 spline-intersection port scope

The **M25 scope study of 11 November 2026**, *M9 and M10 scope: spline
intersections, the port and the proof*, is the source of the historical
measurements and verdict below. Raw evidence retained privately. Its sections
3.2–3.6 give the dependency counts and port units; sections 4.4–4.6 and
Appendix A give the proof experiment; sections 5–6 give the order of work.
This document publishes that scope and separates it from subsequent merged
implementation. It does not declare the full M25 must-show list met.

Evidence cutoff for the implementation references: 12 November 2026, main
`4d51b5c36bdad29632a448585290cfaa627098ba`. The study itself used Catallax
`885db576`, rust-brep main `4c83df59d777ccfcb14f441b1c8d984d1b5d93ab`
(6 November), and NVIDIA-Omniverse/usd-brep SMLib
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`. Its remaining-line counts describe
that historical donor state, not work still outstanding in Catallax today.
Current attribution and pins are in [NOTICE](../NOTICE) and
[THIRD-PARTY.md](THIRD-PARTY.md).

## Four must-show lines

These are the four obligations from [PLAN section 10.1](archive/PLAN-2026-11.md#101-milestones),
preserved verbatim:

1. A written scope of the SMLib intersection subsystem: which SMLib units curve-with-surface and surface-with-surface intersection pull in, directly and through their dependencies, the size of each, and which units serve only SMLib's own tolerance rules and so are not ported.
2. The subsystem cut into port units in dependency order, with the kinds on M8's list first once M8 delivers it.
3. A proof check at scale: on a written case set that includes tangent, touching and near-singular meetings and parts with many spline patches, Catallax proves traced spline intersection curves taken from SMLib's results, and the deterministic work of each proof is published against the 20M application allowance.
4. A written verdict: commit the port order, change it, or narrow what is ported, with the reason. M9's port order and M10's proof route follow the verdict.

## Dependency scope and size

The study built a call graph of the pinned SMLib's 224 C++ source files and
238 headers, then instrumented eleven intersection cases. Physical lines
include comments and blank lines inside definitions; code lines exclude them.
An entered function contributes its whole definition, not just its executed
lines. These counts are neither Rust line counts nor proof-work units.

| Scope at the study pins | Physical lines | Code lines | Already credited in donor, physical / code | Left, physical / code |
| --- | ---: | ---: | ---: | ---: |
| Whole SMLib | 967,172 | 534,277 | 204,203 / 132,995 | — |
| Surface/surface, direct calls only | 99,012 | 64,568 | 27,226 / 18,909 | 71,786 / 45,659 |
| Surface/surface, including reachable virtual overrides | 135,972 | 89,207 | 38,307 / 26,361 | 97,665 / 62,846 |
| Curve/surface, including reachable virtual overrides | 121,983 | 79,643 | 36,611 / 25,078 | 85,372 / 54,565 |
| Entered by the eleven intersection cases | 50,848 | 32,881 | 13,485 / 9,247 | 37,363 / 23,634 |
| Entered by curve/surface | 33,782 | 21,735 | 11,151 / 7,638 | 22,631 / 14,097 |

The virtual-call closure is the planning bound; direct calls alone miss
dynamic dispatch. The surface/surface virtual closure covers 2,253 functions
in 166 files, 14.1% of SMLib's physical lines. The entered set covers 841
functions in 93 files, 5.3%; only 11,899 individual lines executed. Six entered
constructors/conversions were outside the static closure, so it is a bound
within the call-graph model, not an exhaustive compiler-call proof.

The study's family breakdown makes the transitive cost visible:

| Family | Reachable physical lines | Entered physical lines | Entered lines left in donor |
| --- | ---: | ---: | ---: |
| Surface services, cache, point solve and drop | 38,440 | 18,395 | 17,226 |
| Curve services and curve cache | 28,254 | 5,225 | 4,803 |
| NURBS mathematics | 18,839 | 4,688 | 1,087 |
| Topology | 9,830 | 0 | 0 |
| Intersector core | 9,736 | 7,504 | 7,504 |
| Foundation: vectors, boxes, arrays, memory, tolerances | 9,344 | 4,977 | 1,070 |
| Spline objects | 9,036 | 5,490 | 3,046 |
| Solvers and trees | 6,601 | 3,758 | 2,178 |
| Analytic shapes | 5,892 | 811 | 449 |
| **Total** | **135,972** | **50,848** | **37,363** |

The bare-surface cases entered no topology. Trimmed faces, periodic surfaces,
seams and a wider corpus can add dependencies; absence from this entered set
does not authorize pruning them from the general scope. The study found that
the missing layer was mainly objects, caches, search and marching above the
already substantially ported NURBS mathematics, rather than most of SMLib.

### Port units and gates

These are the study's units and historical lines left, in dependency order.
The curve/surface boundary search feeds surface/surface tracing, so Gate A
comes before Gate B. Steps 2 and 3 can follow the foundation independently;
steps 4 and 5 then follow them, with the trim-surface cache needing step 4.
Steps 6–9 follow in order, as recorded by the later
[port foundation PR](https://github.com/BohnBawerick/catallax/pull/388).

| Step | Unit and principal SMLib sources | Physical left | Code left |
| --- | --- | ---: | ---: |
| 1 | Foundation and NURBS helpers: `NL_Util.cpp`, `SmAObject.cpp`, `SmGeomUtility.cpp`, small parts of 30 files | 1,914 | 1,253 |
| 2 | General search: `SmGlobalSolver.cpp`, `SmSolutionArray.cpp` | 1,673 | 1,236 |
| 3 | Spline objects: `SmBSplineCurve.cpp`, `SmBSplineSurface.cpp`, `SmNurbsSrf.cpp` | 2,202 | 1,299 |
| 4 | Curve services and cache: `SmCurve.cpp`, `SmCurveCache.cpp` | 3,313 | 1,916 |
| 5 | Surface cache and point solver: `SmSurfaceCache.cpp`, `SmSurfacePointSolve.cpp`, `SmTrimSrfCache.cpp` | 6,679 | 4,335 |
| 6 | Curve/surface solve, classification and solution combiner: chiefly `SmSurface.cpp` | 6,850 | 4,058 |
| **Gate A** | **Curve/surface differential gate, cumulative** | **22,631** | **14,097** |
| 7 | Marching support: `SmSurface.cpp`, `SmCurve.cpp`, `SmLine.cpp`, `SmBSplineCurve.cpp` and others | 7,120 | 4,805 |
| 8 | Intersector: `SmSurfaceIntersector.cpp`, `SmAdvSurfaceIntersector.cpp` | 7,504 | 4,653 |
| **Gate B** | **Spline surface/plane differential gate, cumulative** | **37,255** | **23,555** |
| 9 | Additional spline/spline path in `SmSurface.cpp` | 108 | 79 |
| **Gate C** | **Spline/spline differential gate, cumulative** | **37,363** | **23,634** |

These are case-driven gates, not acceptance of every path in the static
closure. [M8's intersection list](M8-LATER-CUTS.md#the-intersection-list)
prioritizes fitted spline/plane, then cylinder, cone, sphere and torus, together
with curve/fitted-spline queries. M8 also needs branch coincidence, tangent
germ order and trimmed-face queries. A candidate port alone supplies none of
their Catallax proofs. M8's unsupported analytic pairs remain separate work;
the study establishes no measured completion for them.

### Exclusions and tolerance policy

SMLib's tolerance rules do not become Catallax acceptance rules. Ported
floating-point results remain candidates; Catallax decides completeness,
existence, branch identity and whole-domain deviation under its own rules.
The study's default loft candidate exceeds the existing 1 micrometre cap;
the cap is not relaxed to accept it.

The study does not give a complete per-unit inventory or size for code serving
only SMLib's tolerance policy. That part of obligation 1 remains open.
Nor can all tolerance-related candidate calculations be omitted: the
[merged foundation](https://github.com/BohnBawerick/catallax/pull/388)
includes scalar constants and `SmTol::GetScaledZero` for candidate parity.
That PR replaces `NL_Util` allocators, `SmAObject` attributes/notification and
cache/context bookkeeping with Rust ownership. Those are storage exclusions,
not evidence of a completed tolerance-only exclusion inventory.
[THIRD-PARTY.md](THIRD-PARTY.md) records the actual adapted units and deviations.

## Historical proof-at-scale check

The study measured an exact Python-fractions prototype, **not Catallax
Decider work**. Cells, depths and candidate pieces below are structural
measurements; they cannot be converted into operation units. Every historical
row has **unmeasured work against the unchanged 20M application allowance**.
The reported Python times were indications under load, not benchmarks.

The coverage driver used the named cases `plane`, `tilt`, `rational`, `tight`,
`loop`, `tangent`, `neartangent`, `tangentline`, `cylcut`, `spline` and `curve`.
The wavy and hill inputs are bicubic 6 by 6 nets with nine knot spans each.
Separate `loft` and `lofttight` runs supplied the worked rational patch of
bidegree (2,3), cut by `x + y + z = 5/2`. Nine spans in one patch are not a
measurement of a part with many spline patches.

| Study case | Measured proof structure | Study result or gap | Catallax work / 20M |
| --- | --- | --- | --- |
| Wavy patch, `z = 0.35` | 6 arc, 4 empty cells; depth 1 | Prototype regular-cell cover for two open curves | Unmeasured |
| Hill, `z = 0.3` | 10 arc, 15 empty cells; depth 2 | Prototype regular-cell cover for one loop | Unmeasured |
| Hill, `1e-2` below the peak | 8 arc, 28 empty cells; depth 8 | Small-loop cover | Unmeasured |
| Hill, `1e-4` below the peak | 8 arc, 52 empty cells; depth 14 | Small-loop cover | Unmeasured |
| Hill, `1e-6` below the peak | 12 arc, 76 empty cells; depth 21 | Small-loop cover; SMLib returned no curve | Unmeasured |
| Hill, `1e-8` below the peak | 8 arc, 108 empty cells; depth 28 | Small-loop cover | Unmeasured |
| Hill exactly at peak `225/256` | 64 empty, 4 unresolved cells; depth 16 | Prototype refuses at the touch; SMLib reports success with no curve | Unmeasured |
| Hill, `1e-6` above the peak | 12 empty cells; depth 2 | Prototype proves empty | Unmeasured |
| Rational loft-like patch, tilted plane | 2 arc cells; depth 1; exact ends `(0,0,5/2)` and `(2,0,1/2)` | One regular component; candidate checks below | Unmeasured |
| Cylinder-like patch, tangent plane (`tangentline`) | SMLib returns one tangent curve | No completed prototype proof along the tangent curve | Unmeasured |
| Two wavy spline surfaces (`spline`) | SMLib returns two curves | No spline/spline proof experiment in the study | Unmeasured |
| Parts with many spline patches | No measurement | Scale acceptance remains unshown by this study | Unmeasured |

The loft-like example checks opposite exact signs on both shifted walls of
each cubic parameter-curve piece. A surface-speed bound converts the tube
half-width to a reported 3D bound, taking units as millimetres:

| SMLib candidate | Pieces | Tube half-width `rho` | Reported bound `U` | Against existing limits |
| --- | ---: | ---: | ---: | --- |
| Default setting | 9 | `1/2048` | about `3.806e-3 mm` | Exceeds 1 micrometre cap; refuse or refit |
| Requested tolerance `1e-5` | 21 | `1/131072` | about `5.947e-5 mm` (59 nm) | Numerical bound below 100 nm target, subject to the endpoint gap |

The study checked that the candidate runs monotonically in the free
parameter, but checked its end pieces reaching the patch border by eye.
That missing exact endpoint check prevents treating the prototype tube as a
completed Catallax certificate. Only this worked example checks candidate
tubes; the stress rows measure the cell search. Tangency along a whole curve,
spline/spline proof, many-patch parts and deterministic work were not proved
by those rows. Near-grazing cell counts show increasing subdivision cost;
they do not establish an application-budget pass.

## Verdict and the implementation that follows it

The study's verdict is **keep intersections as the goal, but change the order
of work**:

1. Start M10's plane proof alongside M9, using saved candidates from the pinned
   C++ SMLib. The uncertain proof can be tested before porting the intersector.
2. Direct the port to the intersection chain in the unit order above, Gate A
   then Gate B, then spline/spline. The historical closure is a bounded slice
   of SMLib; the missing object layer is substantial but does not justify
   waiting for the full donor port.
3. Keep the analytic proof order: plane, cylinder, cone, sphere, torus, then
   spline/spline. With an equation on one side, complete search is in two
   patch parameters and must find loops even when the tracer returns none.
   Spline/spline is a separate four-parameter proof problem.
4. Scope free-form fillets separately. Their rails are tangent contacts along
   whole curves, which the crossing-cut certificate does not establish.

The study's proposed first slice required exact ends, a tight candidate within
100 nm, refusal of the coarse and deliberately shifted candidates, recovery
of the missed tiny loop, a named tangent refusal and recorded work within the
unchanged allowance. Failure on ordinary cuts would call for cheaper exact
proof arithmetic before further porting on that proof's account.

Merged evidence shows that work following the verdict, without making the
historical experiment into a completed milestone:

| Follow-on evidence | What it establishes |
| --- | --- |
| [Plane proof PR](https://github.com/BohnBawerick/catallax/pull/367), [committed measurement table](M10-PLANE-SECTIONS.md#verification) | Actual Catallax work including patch admission: tilted loft 96,824; tight candidate 560,119 with `U < 91.871 nm`; default/shifted/short candidates refuse; missed tiny loop 174,443; named tangent refusal 60,998. Exact endpoint checks close the prototype's stated gap. These later counts and bounds are not the study's measurements. |
| [M9 foundation PR](https://github.com/BohnBawerick/catallax/pull/388), [attribution](THIRD-PARTY.md) | The direct Catallax port explicitly follows M25 step 1 and its dependency order, with candidates only and a pinned differential gate. Later units do not retroactively qualify the original proof experiment. |
| [Cylinder/sphere PR](https://github.com/BohnBawerick/catallax/pull/381), [proof record](M10-QUADRIC-SECTIONS.md); [atlas PR](https://github.com/BohnBawerick/catallax/pull/402), [atlas record](M10-ATLAS-SECTIONS.md) | The second slice follows the analytic proof route. Atlas joins, branch identity and third-surface meetings follow study section 4.6; the recorded box/loft Booleans still refuse at later stages. |
| [Fillet scope PR](https://github.com/BohnBawerick/catallax/pull/400), [M27 scope](M27-SPLINE-FILLET-SCOPE.md) | Free-form fillet proof is scoped separately rather than claimed from crossing cuts. |
| [Spline/spline PR](https://github.com/BohnBawerick/catallax/pull/468), [current proof record](M10-SPLINE-SPLINE.md); [R13 decision](R13-SPLINE-TOPOLOGY.md#verdict-for-the-spline-with-spline-stage) | A subsequent 11 November decision selects the 2026 field engine with an exact proof layer. Its wavy pair and loop results carry measured work. This is later proof-route evidence, not a successful spline/spline experiment in the scope study. |

The plane record reports agreement at 1, 2 and 4 workers; the spline/spline
record also checks the default worker count. The requirement remains identical
answers at every thread count. Those bounded observations do not constitute
fresh qualification of every current implementation path.

The [M9 known-miss list](M9-KNOWN-SMLIB-MISSES.md) records R1 exceptions to
candidate parity, including the missed tiny loop and isolated hilltop contact.

## Acceptance still open

The study establishes a written scope, historical unit sizes, structural
proof feasibility and the verdict followed by the port and early proofs.
It does not by itself meet the full must-show list. A complete tolerance-only
exclusion inventory, proof on parts with many spline patches, the full named
tangent/touching/near-singular case coverage and per-proof deterministic work
must still be accounted for against obligation 3. Later measured slices are
partial evidence, not blanket M25 or M10 acceptance.

The 20M application allowance remains a capability requirement. No tolerance,
work charge or allowance is relaxed here. This documentation changes no proven
answer or any of the 477 saved replay answers, and claims no new replay run.
