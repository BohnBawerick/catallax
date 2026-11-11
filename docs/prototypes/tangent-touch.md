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
A subset determinant reserves its coefficient operations through the existing
numeric meter; the meter and all work allowances stay unchanged.
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

**Park.** The round completed 255 requests: 17 cases, five placements and
1, 2 and 4 workers. Every copy had the same answer, including the parameter
enclosures. Worker counts also preserved exact work. No request exhausted 20M.
The example was compiled at release optimisation against the baseline's
numeric and shape libraries. Raw evidence retained privately.

Operation units include patch admission, contact enumeration and the returned
48-bit parameter enclosures. Each number below is identical at 1, 2 and 4
workers. Translation is `(11/3,-17/5,29/7)`, scale is `7/3`, and rotation
uses the exact `3/5,4/5` frame about y. Combined applies all three.

| Case | Outcome | Original | Translated | Scaled | Rotated | Combined |
|---|---|---:|---:|---:|---:|---:|
| bowl | isolated touch | 3,925 | 3,925 | 3,925 | 3,925 | 3,925 |
| saddle | crossing branches | 3,921 | 3,921 | 3,921 | 3,921 | 3,921 |
| flat | degenerate point | 18,240 | 18,240 | 18,240 | 18,240 | 18,240 |
| line | stationary component refusal | 2,462 | 2,462 | 2,462 | 2,462 | 2,462 |
| flush | flush span refusal | 1,630 | 1,630 | 1,630 | 1,630 | 1,630 |
| off-centre | isolated touch | 3,962 | 3,962 | 3,962 | 3,962 | 3,962 |
| several | two isolated touches | 17,839 | 17,839 | 17,839 | 17,839 | 17,839 |
| algebraic | two isolated touches | 24,107 | 24,107 | 24,107 | 24,107 | 24,107 |
| near-miss | no stationary contact | 3,937 | 3,977 | 3,937 | 4,051 | 4,085 |
| tiny-crossing | no stationary contact | 3,925 | 3,959 | 3,925 | 4,039 | 4,067 |
| transverse | no stationary contact | 1,634 | 1,634 | 1,634 | 1,634 | 1,634 |
| boundary | boundary contact refusal | 3,913 | 3,913 | 3,913 | 3,913 | 3,913 |
| hill | isolated touch | 70,566 | 70,566 | 70,566 | 70,566 | 70,566 |
| rational-bowl | isolated touch | 10,699 | 746,611 | 10,699 | 10,699 | 85,653 |
| touch-and-cut | isolated touch beside a regular cut | 9,729 | 9,729 | 9,729 | 9,729 | 9,729 |
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

No repeat or cylinder round was taken. Cylinder is
`prototype.tangent-cylinder-unmeasured`; cone, sphere and torus are respectively
`prototype.tangent-cone-unsupported`, `prototype.tangent-sphere-out-of-scope`
and `prototype.tangent-torus-unsupported`. The experiment has no kernel-loop
request or application consumer, so it cannot qualify Boolean continuation
through those loops. The default plane tangent refusal remains in force.

Reproduce with `cargo run --release -p catallax-intersect --example tangent_touch`.
The example asserts expected successful classifications, copy equality and
worker answer/work equality. Unit checks compare the metered determinant with
the existing fraction-free resultant, verify actual point-on-plane predicates,
and require exhausted work to refuse. No measuring tool, saved answer,
tolerance or work allowance changes.
