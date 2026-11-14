# SMLib intersection port: what it drags in, and whether the proof scales

Outcome: adopted, with the order of work changed.
When: 11 November 2026, published in this repository on 13 November.

## What we were trying to achieve

Cuts and joins on spline faces. The plan was to port SMLib's surface intersection code as a source of candidate curves, and have Catallax prove each curve. Both are now milestone 1.3.

## Why it was proposed

Two doubts stood in the way of committing to that plan, and the owner asked on 11 November for them to be scoped as a task of its own:

1. Does the intersection code pull in most of SMLib, so that the port must nearly finish first?
2. Can Catallax prove traced curves at all, at a cost that fits the allowance?

## What was tried

- A call graph of all 224 SMLib source files, followed from the intersection entry points.
- The C++ SMLib built with coverage and run on eleven intersections, to see which lines really execute.
- The same line ranges compared with what the Rust port already holds.
- The missing proof worked by hand and by an exact-arithmetic script, on a loft-like patch cut by a plane, using a curve traced by the C++ SMLib.

## Result

| Question | Number |
|---|---:|
| Lines the intersector can reach, worst case | 135,972 (14.1% of SMLib) |
| Lines that ran for the eleven intersections | 50,848 (5.3%) |
| Of those, already ported | 13,485 |
| Still to port before a spline face can be cut by a plane | 37,255 |
| Lines the six open port changes added toward it | 0 |

The proof goes through cheaply for a clean crossing. It gets slower near a grazing contact and refuses at an exact tangent contact.

One finding nobody expected: for a plane one millionth of a unit below a hilltop, SMLib returned no curve. The exact check found the small closed loop.

## Why it was adopted, and what changed

Intersections stayed the goal. Three things changed:

- The proof started first, on the plane cut, with curves saved from the C++ SMLib. It did not wait for the port.
- The port was pointed at the intersection chain.
- Catallax keeps its own complete search whatever the tracer says, because the tracer can miss a piece. For a spline face cut by a plane, cylinder, cone, sphere or torus that search needs no tracer at all. The tracer earns its place only for spline against spline. SMLib's own misses are listed as known misses, so a comparison with SMLib does not reward them.

## Sources

- [M25-INTERSECTION-PORT-SCOPE.md](../M25-INTERSECTION-PORT-SCOPE.md): https://github.com/BohnBawerick/catallax/pull/538
- [M9-KNOWN-SMLIB-MISSES.md](../M9-KNOWN-SMLIB-MISSES.md): https://github.com/BohnBawerick/catallax/pull/563
- [archived plan](../archive/PLAN-2026-11.md) section 12.6, row R1
- The scope study itself, "M9 and M10 scope", 11 November 2026. Its raw evidence is kept outside this repository.
