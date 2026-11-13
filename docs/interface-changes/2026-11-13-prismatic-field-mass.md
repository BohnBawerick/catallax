# Exact field CSG moments and one reported mass

Pieces: field mass and application assembly queries. Additive opt-in APIs:
`mass::prismatic_properties`, `assembly::prismatic_uniform_mass` in the field,
and `assembly::resolve_mass` in the application. Existing mass queries,
recipe formats, manifests and export routes are unchanged.

## Covered geometry

The new field query accepts bounded coordinate-axis boxes and finite right
circular cylinders parallel to coordinate axes. A planar extrusion qualifies
only when its complete convex-plane definition proves it is an axis box.
A halfspace cell qualifies only when every defining constraint is a rational
axis plane, with a finite lower and upper bound on each coordinate. Contradictory
or zero-width cell bounds contribute zero material, including in subtraction.
An enclosing box alone cannot replace a curved or concave leaf.

Union, intersection and subtraction use exact indicator algebra:
`A union B = A + B - A*B`, `A intersect B = A*B`,
`A subtract B = A - A*B`. The resulting signed intersection terms integrate
the original field material. Bounds remove a term only when its intersection
has zero volume. This is a mass proof, not a contact or clearance classifier.

Each remaining intersection must contain parallel cylinders. It is then a
Cartesian product of a finite axial interval and a rectangle intersected with
disks. Disjoint components may have different axes. A nonzero candidate with
nonparallel cylinder axes refuses `field.mass-prismatic-coverage` with its
enclosing box. Curved profiles, non-axis frames, cones, lattices and other
uncovered leaves also refuse. No approximate geometry substitutes for them.

Equal-radius disks with collinear centres retain only their extreme centres
in an intersection. Squared distance is convex in the centre, so those two
constraints imply all intermediate disks. The complete intersection bounds,
including each finite cylinder's axial limits, remain in the term.

## Boundary integration

Every circle/circle and circle/rectangle carrier intersection supplies an exact
vertex. Exact predicates order and deduplicate vertices. Between consecutive
vertices, an open boundary cell crosses no other carrier. One exact point
therefore proves that entire cell inside or outside all the other constraints.
An unresolved predicate or a point on another carrier refuses; it never becomes
a sampled membership guess.

Accepted circle arcs run counterclockwise. The two vertical rectangle sides
have their usual boundary orientations; horizontal sides have `dy = 0`.
Green's identity
`integral x^p*y^q dA = integral_boundary x^(p+1)*y^q dy/(p+1)`
gives all six section moments through degree two. Straight segments integrate
algebraically. Circular terms use algebraic endpoint expressions and certified
arc angles from checked arctangent and pi enclosures. Axial antiderivatives
then give all ten spatial coordinate moments. The existing moment-to-centroid
and inertia calculation supplies the requested world properties.

Box-only CSG remains rational throughout, including centroid and inertia. It
does not lose exact rational answers through a dyadic enclosure, and can meet
zero requested width. A regression reproduces and closes that unnecessary
precision refusal.

The existing cell limit caps processed CSG intersections and boundary cells,
including refinement rounds. Work is charged to the caller's active arithmetic
allowance. No query creates a fresh allowance. Numeric, cell, precision and
coverage failures retain typed refusals and the counts or enclosure reached.
An incomplete sum is never reported as a volume enclosure.

## Complete authored P1

The full public bracket now has an independent field mass and centre of gravity
inside the unchanged 20,000,000-unit application allowance. The original field
mass route still refuses `field.print-coverage`, and B-rep construction keeps
its M2 refusal. The new route evaluates every authored cut, including all eight
overlapping slot cylinders and the slot's rectangular middle.

An independent component formula checks the result:

```text
V = 19080 + 275*pi                         mm3
Mx = 30*V                                 mm4
My = 257700 + (16255/2)*pi                 mm4
Mz = 153690 + 3590*pi                      mm4
centre = (Mx/V, My/V, Mz/V)                mm
```

At explicit density `0.00785` g/mm3 the answer is about 156.559913140937 g,
with centre near `(30,14.201472877614,8.271602017314)` mm. Query work is
2,795,348 at 1, 2 and 4 threads. Exact enclosures and work compare equal before
decimal report rounding. The query also meets the unchanged inertia accuracy;
no unrequested-property tolerance is relaxed to admit the mass result.

Analytic controls check full cylinder and clipped-quarter-disk moments,
nonzero mixed inertia, equal- and unequal-radius lens intersections, unions
and differences, a capsule with redundant intermediate drills, different-axis
disjoint components, empty material, exact rational box CSG and exhausted
resource controls. Nonparallel overlapping cylinders retain a typed refusal.

## Complete authored P4

The tray's compiled shell contains the original block minus a planar cavity.
The new query proves the cavity's complete halfspace definition is an axis box,
then integrates the authored radius-23 clipping cylinder. It does not rebuild
the shell from its display or an envelope.

A separate circular-segment formula checks the volume and first z moment.
Let `A` be the area of the centred 40 by 30 rectangle inside the radius-23 disk.
The inner 36 by 26 rectangle lies wholly inside that disk.

```text
A = 2*(15*sqrt(304) + 20*sqrt(129)
       + 529*(atan(20/sqrt(129)) + atan(15/sqrt(304)) - pi/2))
V = 12*A - 9360                            mm3
Mz = 72*A - 65520                          mm4
centre = (0,0,Mz/V)                        mm
```

At density `0.00785` g/mm3 this gives mass
`[37.92534104783266,37.92534104783267]` g and centre z
`[4.062614653686839,4.062614653686840]` mm. Complete compilation, placement and
measurement use 40,411 work units at 1, 2 and 4 workers, with identical exact
enclosures and work. The default inertia accuracy also passes. The previous
integrator exhausted its unchanged 2,000,000-cell limit with volume width
`158617329/262144` mm3, about 605.077091, and maximum centroid width
`66519527507/51035508096` mm, about 1.303397. That route remains unchanged.

## Reporting one value

Application `resolve_mass` consumes already completed independent query
results. It performs no build and creates no work allowance.

- Two proven results must pass the existing bounded hull agreement check.
  Disagreement or an excessively wide hull refuses rather than selecting one.
- One proven result is checked at the requested accuracy and returned with
  `MassSource::Brep` or `MassSource::Field`. `ResolvedMass::unavailable` keeps
  the other route's original refusal, including its undecided numeric bound.
- Two refusals return `assembly.mass-unresolved` with both reasons and the
  first undecided decision retained. Callers retain both complete query records.

Only `MassSource::Agreed` claims independent agreement. A field-only lattice
mass is valid even though its B-rep region envelope cannot supply material
mass. A B-rep-only proven value remains available when field integration
refuses. [The M17 record](../M17-RELATION-MASS.md) lists each target's source.

The full gate now requires one reported mass/centroid on P1-P4 and P6, and pins
authored P5's named refusal. It compares unrounded results, refusals, undecided
numeric bounds, selected sources and operation work at 1, 2 and 4 threads.
P5 is not approximated by a smaller fillet or a surrogate body. Raw evidence
retained privately.

## Remaining coverage

Independent field integration of P3 remains open; its B-rep value is the
reported authority. General oblique-frame normalisation, nonparallel
cylinder intersections, curved prismatic profiles and heterogeneous density
are further work. These limits do not change the caller's accuracy, material
model, recipe or field definition.
