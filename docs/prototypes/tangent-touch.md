# Exact plane tangent contacts

## Question and bounded round

R2 asks whether exact elimination can enumerate nonregular contacts of a
plane with a retained rational patch inside the unchanged 20,000,000 work
allowance. The baseline is `eee8fd123ce8a97e1115013d436f9db7e51d3043`.
M10 refuses its hill tangent as `intersect.patch-plane-tangent`.
This file precedes the experiment's code. Run one measured round, with at most
one repeat if the measurement is wrong. No wall-clock target applies.

## Candidates and cases

Compare the existing complete cell search with exact elimination of the two
derivatives of the homogeneous plane equation F. Isolate roots of both
projection resultants, check every coordinate pair against F, Fu and Fv,
and classify the surviving points by the exact Hessian determinant.
The subset determinant and the equation and derivative producers reserve their
coefficient operations through the existing numeric meter; the meter and all work allowances stay unchanged.
A positive determinant proves an isolated touch; a negative determinant
proves two crossing branches. A zero determinant refuses as degenerate.
A zero projection, a flush span, a nonsmooth knot join or exhausted work
must refuse, never claim an empty contact set.

The cases include the M10 hill tangent, a bowl, a saddle, a flat higher-order
touch, a contact line, a flush face, an off-centre touch, several touches,
an algebraic parameter touch, a touch beside a regular cut, a boundary
contact, a creased patch, and separated and transverse cuts. Run
translated, uniformly scaled and rationally rotated copies, including
combined copies, with 1, 2 and 4 independent request workers. Publish each
case's operation work, including patch construction, against 20M. Record
cylinder, cone, sphere and torus coverage explicitly.

The polynomial graph fixtures use height F/16 so the existing patch validator
can certify their whole-domain chart. Keep the unscaled off-centre graph as
`steep-off-centre`, and retain the creased graph, to measure admission refusals
separately from elimination. The M10 hill fixture keeps its original height.

## Decision rule

Promote only if every required plane case has a complete typed classification
or a proved named degeneracy within 20M, with identical classifications on
all copies and identical answers and work at 1, 2 and 4 workers. Otherwise
park. Finite contact enumeration alone does not certify a Boolean: the
complete regular cut and the isolated points must still share a coverage
certificate. Promotion needs a separate production implementation.

The dispatch requires branch `fm/cx-r2-tangent-touch` and a reviewable PR,
overriding section 8.4's prototype branch name. Keep the experiment in an
opt-in example, outside the default proof path. The application loop cannot
consume these contact reports. Run the existing plane gate as the baseline
and leave saved replay inputs, expected answers and measuring tools unchanged.

## Result

**Park.** The corrected round completed 255 requests: 17 cases, five placements and
1, 2 and 4 workers. Every copy had the same answer, including the parameter
enclosures. Worker counts also preserved exact work. No request exhausted 20M.
The example ran through Cargo at release optimisation with the baseline's
unchanged numeric and shape code. Raw evidence retained privately.

The first round at `e2c9285a` omitted eager rational work in plane-equation
assembly, first and second derivatives, and determinant sign changes. Its
70,566-unit hill count and 746,611-unit maximum are withdrawn. The single
allowed repeat corrects that measurement: each coefficient operation now
reserves work through `Decider::reserve_exact_arithmetic` before running.
All 255 answers and parameter enclosures match the first round exactly.

Operation units include patch admission, contact enumeration and the returned
48-bit parameter enclosures. Each number below is identical at 1, 2 and 4
workers. Translation is `(11/3,-17/5,29/7)`, scale is `7/3`, and rotation
uses the exact `3/5,4/5` frame about y. Combined applies all three.

| Case | Outcome | Original | Translated | Scaled | Rotated | Combined |
|---|---|---:|---:|---:|---:|---:|
| bowl | isolated touch | 4,010 | 4,055 | 4,010 | 4,090 | 4,135 |
| saddle | crossing branches | 4,006 | 4,051 | 4,006 | 4,086 | 4,131 |
| flat | degenerate point | 18,325 | 18,370 | 18,325 | 18,405 | 18,450 |
| line | stationary component refusal | 2,517 | 2,562 | 2,517 | 2,572 | 2,617 |
| flush | flush span refusal | 1,665 | 1,710 | 1,665 | 1,695 | 1,740 |
| off-centre | isolated touch | 4,102 | 4,142 | 4,102 | 4,217 | 4,252 |
| several | two isolated touches | 18,174 | 18,214 | 18,174 | 18,354 | 18,389 |
| algebraic | two isolated touches | 24,247 | 24,287 | 24,247 | 24,377 | 24,412 |
| near-miss | no stationary contact | 4,031 | 4,114 | 4,031 | 4,257 | 4,332 |
| tiny-crossing | no stationary contact | 4,019 | 4,096 | 4,019 | 4,245 | 4,314 |
| transverse | no stationary contact | 1,709 | 1,754 | 1,709 | 1,749 | 1,794 |
| boundary | boundary contact refusal | 4,023 | 4,063 | 4,023 | 4,113 | 4,148 |
| hill | isolated touch | 71,036 | 71,056 | 71,036 | 71,221 | 71,231 |
| rational-bowl | isolated touch | 10,854 | 746,891 | 10,854 | 10,954 | 86,028 |
| touch-and-cut | isolated touch beside a regular cut | 9,929 | 9,974 | 9,929 | 10,059 | 10,104 |
| crease | chart admission refusal | 2,936 | 2,936 | 2,936 | 2,936 | 2,936 |
| steep-off-centre | chart admission refusal | 155,393 | 155,393 | 155,393 | 279,485 | 279,485 |

The baseline plane gate passed all four tests, including its 1/2/4 worker
comparison. It refused the original hill tangent at 60,998 units. Elimination
finds its exact parameters `(3/2,3/2)`; the algebraic case finds
`(+/-sqrt(1/2),0)`. The off-centre point is `(1/3,2/5)`.

The finite interior classification works. Full R2 does not meet promotion:
`prototype.tangent-stationary-component` does not resolve positive-dimensional
critical sets, and `prototype.tangent-boundary-contact` leaves an ordinary
boundary touch undecided. A zero Hessian returns the named local refusal
`prototype.tangent-degenerate-contact`; a node returns
`prototype.tangent-self-touch`. Flush spans refuse as
`prototype.tangent-flush-span`. The two admission cases retain
`shapes.rational-patch-unproven`. These are reported gaps, not decided empty
cuts. A report with no stationary contacts says nothing about regular cuts.

The one allowed repeat was used for the metering correction. No cylinder
round was taken. Cylinder is
`prototype.tangent-cylinder-unmeasured`; cone, sphere and torus are respectively
`prototype.tangent-cone-unsupported`, `prototype.tangent-sphere-out-of-scope`
and `prototype.tangent-torus-unsupported`. The experiment has no kernel-loop
request or application consumer, so it cannot qualify Boolean continuation
through those loops. The default plane tangent refusal remains in force.

Reproduce with `cargo run --release -p catallax-intersect --example tangent_touch`.
The example asserts expected successful classifications, copy equality and
worker answer/work equality. Unit checks compare the metered determinant with
the existing fraction-free resultant, verify actual point-on-plane predicates,
and require zero-allowance equation assembly, derivatives and sign changes
to refuse before arithmetic. The original zero-allowance equation test
failed because construction succeeded without reserving work; all seven
example tests now pass. No measuring tool, saved answer,
tolerance or work allowance changes.
