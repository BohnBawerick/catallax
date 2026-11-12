# Exact field moments of coaxial polygonal revolutions

Piece: field mass. Two opt-in entries are added:
`mass::revolved_properties` and `assembly::revolved_uniform_mass`.
The existing `mass::properties` and `assembly::uniform_mass` routes are
unchanged. No recipe, manifest, STEP or replay format changes.

## Coverage and proof

Every live leaf must retain a complete, full revolution with a polygonal
meridian. Axes must be coaxial. Their unit directions, meridian vertex heights
and radii must be rational. Loops can have holes; Boolean union, intersection
and subtraction use the original field instruction tape. Rigid placements keep
their world coordinates. Curved profiles, partial turns, noncoaxial axes and
uncovered irrational frames refuse `field.mass-revolution-coverage`.

The algorithm splits the axis at every profile vertex height and every
intersection of active affine radial boundaries. Within each open axial slab,
the radial boundaries have a fixed order. Consecutive boundaries delimit a
connected cell that meets no leaf boundary. An exact field membership query
therefore classifies that entire cell. Unresolved membership refuses; it is
never inferred from a display sample. Caps, rails and the axis have zero volume.
Outside the finite meridian bounds every leaf is outside, so the finite cells
cover all material of every Boolean combination.

For a material cell with inner radius `a(h)` and outer radius `b(h)`, both affine:

- `V/pi = integral (b^2 - a^2) dh`.
- `M/pi = integral h (b^2 - a^2) dh`.
- `J_axis/pi = integral (b^4 - a^4) dh / 2`.
- `J_perpendicular/pi = J_axis/(2*pi) + integral h^2 (b^2-a^2) dh - M^2/(pi*V)`.

Polynomial antiderivatives give exact rational coefficients. The axial centre
`M/V` is rational. The world inertia tensor is
`J_perpendicular * identity + (J_axis-J_perpendicular) * u*u^T`, where `u` is
the exact unit axis. Certified pi enclosures supply volume and inertia bounds.
Precision increases only under the caller's existing numeric policy, until
every requested width is met or the query refuses. The rational property
proposal lies inside those bounds.

The existing cell cap limits meridian edges and processed cells. Vertex,
boundary-pair, tape and integration work is charged to the active operation.
No new allowance replaces the caller's request. Exhaustion retains typed
numeric or `field.mass-budget` refusals with counts.

## Authored P2

The complete collar and both outer chamfers now have an independent field mass
and centre of gravity at default accuracy. The unchanged original field mass
route still refuses `field.print-coverage`; the new query does not relabel that
route as covered. The application probe explicitly tries the opt-in meridian
query and uses the original integrator only for uncovered definitions.

For the authored dimensions, the annulus contributes `1296*pi` mm3. Each
one-millimetre chamfer removes `35*pi/3` mm3, giving `3818*pi/3` mm3 and centre
`(0,6,0)` mm. At the explicit test density `0.00785` g/mm3, the mass is
`299713*pi/30000` g, about 31.385872 g. B-rep and field enclosures pass the
existing hull agreement check, so the consumer receives one bounded value.
The new field measurement spends 389 work units at 1, 2 and 4 threads, with
identical exact enclosures, reported values and work. Compilation and placement
are separate from that query count.

Controls check annulus volume and all inertia components against analytic
values, Boolean meridians whose radial boundaries cross between profile
vertices, a nested ring cavity, quarter-turn and oblique placements, empty
material, uncovered partial/noncoaxial definitions, invalid density and
exhausted work/cell limits. The complete P2 test also pins the old route's
coverage refusal.

P1, P3, P4 and P6 do not fit this route. P5 retains its compilation refusal.
[The M17 record](../M17-RELATION-MASS.md) states their independent mass outcomes
and remaining limits. This change does not complete M17. Raw evidence retained
privately.
