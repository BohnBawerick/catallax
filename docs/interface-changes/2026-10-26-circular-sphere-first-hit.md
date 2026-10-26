## 2026-10-26: Exact circular first-hit profiles on stock spheres

Pieces: recipe features, field; application, measurement, display and STEP consumers.

Change: `AuthoredCap::compile` also admits one full rational circular sketch
against one stock sphere when its centre lies on the sphere/profile axis and
the raw sweep direction is parallel to that axis. Both `extrude_to` on the
stock wall and `extrude_next` select the same certified positive entry or exit.
The existing polygon and annular-torus routes retain their behavior.

Why: a complete disk has an exact cylindrical side and spherical end. A polygon
approximation would change its definition, volume, rim and authored origin.

Let `R` be the sphere radius, `r` the disk radius, `e` its unit normal oriented
along the ray, `s > 0` the raw ray speed and `z0` the signed start height relative
to the sphere centre. Compilation proves coaxiality and `0 < r < R` exactly.
For every `0 <= rho² <= r²`, the two ray roots are
`(-z0 ± sqrt(R²-rho²))/s`, with discriminant at least
`4 s² (R²-r²) > 0`. Entry requires `z0 < -R`; exit requires
`z0² + r² < R²`. These strict conditions include the axis and the closed rim,
exclude every zero-hit or changing-branch boundary, and give a uniform positive
lower bound: `(-z0-R)/s` for entry and
`(-z0+sqrt(R²-r²))/s` for exit. Every selected root is bounded above by
`(-z0+R)/s`. Reversing or scaling the authored ray preserves this definition.

The B-rep directly builds the exact plane, cylinder and selected spherical cap
with two shared ring edges. The end is outward for exit and inward to the
sphere for entry. Its spherical chart aligns with the circular sweep axis,
making the rim an exact constant-latitude trim. The pole belongs to the cap's
interior and needs no extra topological vertex; no geometry is perturbed.
The independent field compiler consumes a finite cylinder/plane column and
the original spherical implicit definition. Exit intersects the column and
sphere. Entry intersects the column with a finite midpoint halfspace and
subtracts the sphere. Every artificial bounding face is excluded by the same
certificate. No B-rep topology or containment supplies field geometry.

`start`, `side.segment` and `end` belong to the extrusion. The circle's authored
id names its one cylindrical wall; cap and start edges retain the corresponding
physical origin names through radius and start-height edits.

Migration: existing polygon callers need no changes. The new
`CurvedCapped::circular_radius()` returns `Some(r)` for a disk, whose chart origin
is its circle centre and whose exact domain is `u²+v² <= r²`. Its linear
`side_constraints()` is empty and `side_ids()` contains the one circle id.
Field membership and origin queries consume this radial constraint explicitly.
Construction cells and the public recipe formats are unchanged.

Remaining coverage: silhouettes and tangencies, zero-hit/branch-changing starts,
offset axes, oblique directions, annuli on spheres, general curved profiles and
other target graphs retain coverage refusals. No cone-apex, native-display,
print-route, exporter, tolerance, allowance or resource-cap change is included.
Native field display and requested field mass retain their independent consumer
outcomes; construction and sampled agreement do not certify those outputs.
