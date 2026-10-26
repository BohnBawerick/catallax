# A planar opening below a spherical cap

The ordinary recipe evaluator and independent field compiler accept one convex
coaxial cell with a cylindrical wall, one spherical cap and one perpendicular
planar base selected open. The base must lie between the sphere centre's plane
and the cylindrical/spherical join, including the centre plane and excluding the
join. Either axial direction is accepted. The input must have complete circular
trims and outward material senses. This extends the distance-defined shell in
[RECIPE 4.7](RECIPE.md#47-shell) within the existing circular cavity representation.

The public [fixture](../crates/catallax-features/tests/fixtures/open-spherical-cap/base-open.json)
intersects a radius-3 cylinder at `0 <= z <= 6` with a radius-5 sphere centred at
the origin, then shells by `1/2` with `stock:start` open. The actual source join
is `z=4`, and its top pole is `z=5`. The inner radius is `5/2`, its spherical radius
is `9/2`, and its join is `z=sqrt(14)`. The original planar base becomes an annulus
with radii `5/2` and `3`. There is no lip beyond the original solid.

## Whole-domain finite-face proof

Orient the coordinate `u` from the open base towards the sphere's pole, with the
sphere centre at `u=0`. Let the cylinder radius be `a`, sphere radius `R`, base
height `b`, and join height `h=sqrt(R²-a²)`. The accepted hypotheses are
`R>a>0` and `0<=b<h`. Every source point has radial coordinate `0<=r<=a`,
`u>=b`, and `q=sqrt(r²+u²)<=R`. The closed kept faces are the cylinder at
`b<=u<=h` and the sphere at `h<=u<=R`.

Let `dc=a-r` and `ds=R-q`. Each is a lower bound on distance to its finite face.
The following proves that the smaller bound has a witness on its finite face:

1. If `dc<=ds`, then `q-r<=R-a`. If `u>h`, squaring gives
   `sqrt(r²+h²)>=r+R-a` for `r<=a`, and hence `q-r>R-a`, a contradiction.
   The cylindrical projection has `b<=u<=h`. On the axis any azimuth supplies
   a witness. At `u=b` that witness lies on the closed rim circle.
2. If `ds<=dc`, then `r<=q-(R-a)<=a*q/R`, since `q<=R` and `R>a`.
   As `u>=0`, this implies `R*u>=h*q`. The radial spherical projection is
   therefore on the finite cap. `q=0` cannot minimise `ds`, because `a<R`.
3. The non-strict inequalities include exact ties and the join. Both source
   faces include their boundary circles.

Thus `d(p,K)=min(dc,ds)` throughout the source. The cavity is exactly the source
base halfspace intersected with `r<a-t` and `q<R-t`, regularised on subtraction.
Only after checking the hypotheses and physical section inventory does
`CircularSkin::cavity` use these inward supports. The B-rep recogniser independently
checks complete source loops, traversal and material senses; the field extracts
supports from the recipe instruction tree without consulting B-rep topology.

If the inward cylinder is consumed, the source is retained. A raised base can
consume the inner cylindrical band, including an exact join/base tie, while
leaving a spherical cavity; an axially consumed cavity also retains the source.
The existing consumers return empty when every source face is selected open.

## Defect and control ledger

| Case | Evidence and required behavior |
|---|---|
| Base-open fixture | Reproduced on `866c9a30c4730ebe906cba352f61bf350b251652`: the real verification recipe and field entries returned `distance.cylinder-domain`; the behavioral build test failed at that refusal. The same fixture now builds in both kernels. |
| Finite distances | An independent rational oracle clamps the cylindrical projection to `[0,4]` and replaces an out-of-cap spherical projection with distance to the join circle. It checks 375 points, plus exact wall, rim, join and pole controls. This sampling supplements the whole-domain proof above. |
| Thickness and material | Probes straddle the inner cylinder and sphere, include the open hole and exterior below the rim, and distinguish a full-sphere support point below the inner join from an exposed inner face. B-rep material classification and independent field classification are exercised separately. |
| Origins and finite sections | The source has three faces and two complete circles. The result has five faces, four complete circles and no invented pole vertices. Inner faces carry `inner.of`, kept faces retain source names, and the base annulus retains `stock:start`. Inner and outer material senses are checked. |
| Independent properties | Exact disk integrals enclose volume, area, centroid and centroidal inertia under the default property accuracy. For the fixture, volume is `pi*(112*sqrt(14)-241)/12` and area is `pi*(309/4-4*sqrt(14))`. |
| Poses and scale | Rational non-cardinal rotation, translation, an axis reversal, scales `1/10` and `10`, and a decreasing-axis cap preserve the finite-face material definition. |
| Topology transitions | Exact inner join/base tie, disappearance of the inner cylindrical band, radial and axial cavity consumption, and selection of every face retain regularised semantics. These are controls, not separate reproduced pre-change defects. |
| Application and output | The real application entry and an appended shell from a saved source prefix use the existing 20M work allowance and default mesh controls. Display triangles form a closed oriented mesh; fresh and prefix builds retain identical names and mesh. STEP is imported without repair, with independent volume and hole/material probes. |
| Unsafe extension control | A base at `z=-4` stays refused. At `(0,0,-4)` the full sphere is distance `1` away but the finite kept cap is distance `sqrt(73)` away; the actual nearest kept cylinder is distance `3` away. A support minimum would give wrong material at thickness `3/2`. Opening the spherical face also stays refused. |

## Remaining coverage

This is not a general open spherical-shell envelope. Bases below the centre
plane, several spherical bands, partial circles, noncoaxial supports, reversed
source material and open curved faces do not acquire this proof. Cases needing
finite rim pipes continue to report their existing coverage refusal. No torus
envelope representation, dispatcher, native display implementation, STEP writer,
accuracy target, deterministic allowance or independent vertex limit changes.
Native field display and sampled background comparison remain separate consumers;
their outcomes do not replace the whole-domain proof or independent properties.

The tests are [finite shell controls](../crates/catallax-features/tests/open_spherical_cap.rs),
[application controls](../crates/catallax-app/tests/open_spherical_cap.rs), and
[independent STEP controls](../tools/step-regressions/test_open_spherical_cap.py).
Public closed circular-shell fixtures are replayed against the preserved baseline
with strict manifest and STEP-byte comparison. Final release preservation remains
subject to the existing release policy and commit-matched verification.
