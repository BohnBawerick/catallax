# M21: imported spheres and tori cut about two axes

## What main already supplied

The starting commit `015f64b8` already built complete spheres and ring tori,
planar parts, single cylindrical or conical walls between disks, pointed cones,
and one-axis bodies of revolution. The latter include hemispheres, domes,
lenses, spherical zones, torus halves and caps, and rounded cylinders.
`crates/catallax-compare/tests/import_revolution.rs` checks their face names,
membership, volume and requested agreement. The
[one-axis interface note](interface-changes/2026-11-12-field-step-revolved-spheres-and-tori.md)
records that boundary. Those successful routes still run first.

Flat bilinear spline patches and saddle roofs also already build. Their
tests are in `crates/catallax-compare/tests/import_spline.rs`. This change
does not rebuild that work.

## Added analytic definitions

After the revolution reader declines a domain, the field can retain:

- One outward sphere closed by planar disks about different axes. Every cut
  removes a minor cap, and every disk must lie strictly inside all the other
  cuts. The admitted sphere has radius 5 and two perpendicular cuts at height
  4, each with a radius-3 rim. Its volume is `472 pi / 3`.
- One outward ring torus closed by two meridian disks, bounding a convex sector
  strictly between 0 and pi. The admitted torus has major radius 4 and minor
  radius 1. Tests cover a right-angle sector, an acute sector with radial
  directions `(1,0)` and `(3/5,4/5)`, and an obtuse sector ending at
  `(-3/5,4/5)`. Their volumes are `2 pi^2`,
  `4 pi atan(4/3)` and `4 pi (pi - atan(4/3))`.

The field retains an intersection of the exact analytic supports. It keeps
every admitted face name and no imported topology for subsequent membership
queries. Every boundary is a shared whole circle, traversed in opposite
directions by its two faces. The reader checks disk orientation, incidence,
radius and closure before publishing a field.

For sphere cuts, let `g` be the distance from one disk centre to another
outward plane, `rho` its rim radius and `c` the two unit normals' dot product.
Admission requires `g > 0` and `g^2 - rho^2 (1-c^2) > 0`. Both values appear
in a refusal if separation fails. This proves that removed caps do not overlap
and that their volumes can be subtracted independently.

For torus sectors, the two disk centres must lie exactly on the major circle,
their radii must equal the tube radius, and their planes must contain the
torus axis. Each radial centre direction must lie strictly inside the other
plane. These signs select the convex sector and exclude the other intersection
circle of each complete plane. Volume uses a checked arctangent enclosure
under the existing accuracy, precision, deadline and work controls.

## Remaining named refusals

| Case | Status and refusal | What resuming requires |
| --- | --- | --- |
| Crossing or touching spherical caps, split circles, mixed curved carriers | `field.import-domain`, with face or arc counts, or the failed separation gap and squared margin | Finite arc domains and their coverage proof |
| Torus sectors with coplanar disks or a reflex angle | `field.import-domain`, with cross-plane signs. A half torus of radii 4 and 1 has signs `0, 0` | Distinguishing multiple disk domains on the same support, or a union of angular cells |
| General spline faces | Parked behind M10. `field.import-spline` remains the named refusal. Flat equal-weight degree-1 patches and six-face saddle roofs are the existing exceptions | M10's finite spline domains, membership and complete boundary coverage in the field |
| M16 solver queries on imported parts | Parked by the release decision; no new solver-query claim | A separate assignment reopening that scope |

The new tests exercise the real STEP reader, shared admission, field compilation,
face selection, exact membership and volume, and requested two-kernel agreement.
They repeat placed copies and compare answers and numeric work with 1, 2 and 4
field workers. Reversed curved faces and reversed disks must refuse. A clipped
torus volume must refuse when its work allowance or deadline is spent.

No tolerance or work allowance changes. The 477 saved answers remain the
preservation gate. Raw evidence retained privately.
