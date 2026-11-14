# SMLib intersection port: what it drags in, and whether the proof scales

Outcome: adopted, with the order of work changed.
When: 11 November 2026, published in this repository on 12 November.

## What we were trying to achieve

Cuts and joins on spline faces. The plan was to port SMLib's surface intersection code as a source of candidate curves, and have Catallax prove each curve. Both are now milestone 1.3.

## Why it was proposed

Two doubts stood in the way of committing to that plan, and the owner asked on 11 November for them to be scoped as a task of its own:

1. Does the intersection code pull in most of SMLib, so that the port must nearly finish first?
2. Can Catallax prove traced curves at all, at a cost that fits the allowance?

## What was tried

- A call graph of all 224 SMLib source files, followed from the intersection entry points.
- The C++ SMLib built with coverage and run on eleven intersections, to see which functions they enter and which lines execute.
- The same line ranges compared with what the donor Rust port (rust-brep at the study pin) already held.
- The missing proof worked by hand and by an exact-arithmetic Python prototype, on a loft-like patch cut by a plane, using a curve traced by the C++ SMLib. The same prototype ran a cell search on wavy and hill patches cut by planes.

## Result

All line counts are physical lines of C++ source, comments and blank lines inside definitions included. An entered function counts its whole definition.

| Question | Number |
|---|---:|
| Physical lines in the surface-with-surface call graph, virtual calls included. This is the modelled planning bound. | 135,972 (14.1% of SMLib) |
| Physical lines of the 841 functions the eleven intersections entered | 50,848 (5.3% of SMLib) |
| Individual lines that actually executed in those runs | 11,899 |
| Of the 50,848 entered lines, already in the donor port | 13,485 |
| Entered lines still to port before a spline face can be cut by a plane (gate B) | 37,255 |
| Entered lines still to port for spline against spline (gate C) | 37,363 |

The prototype's cell search covered a clean crossing with few cells: 2 to 10 arc cells at depth 1 or 2. Near a touch it needed many more and deeper cells: depth 21 for a plane 1e-6 below a hilltop, depth 28 at 1e-8. Exactly at the touch it refused. The study did not measure Catallax work units, so none of these runs is counted against the 20M allowance.

One finding nobody expected: for a plane 1e-6 below a hilltop, SMLib returned no curve. The exact search found the small closed loop.

## Why it was adopted, and what changed

Intersections stayed the goal. Three things changed:

- The proof started first, on the plane cut, with curves saved from the C++ SMLib. It did not wait for the port.
- The port was pointed at the intersection chain.
- Catallax keeps its own complete search whatever the tracer says, because the tracer can miss a piece. For a spline face cut by a plane, cylinder, cone, sphere or torus, the complete search runs in the two patch parameters and must find loops even when the tracer returns none. Spline against spline is a separate proof problem in four parameters. SMLib's own misses are listed as known misses, so a comparison with SMLib does not reward them.

## Sources

- [M25-INTERSECTION-PORT-SCOPE.md](../M25-INTERSECTION-PORT-SCOPE.md): https://github.com/BohnBawerick/catallax/pull/538
- [M9-KNOWN-SMLIB-MISSES.md](../M9-KNOWN-SMLIB-MISSES.md): https://github.com/BohnBawerick/catallax/pull/563
- [archived plan](../archive/PLAN-2026-11.md) section 12.6, row R1
- The scope study itself, "M9 and M10 scope", 11 November 2026. Its raw evidence is kept outside this repository.
