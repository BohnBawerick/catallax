## 2026-11-13: the field reads two-axis sphere cuts and torus sectors

Pieces: field, compare

`Plan::compile_with_files` now retains spheres with disjoint minor caps removed
by nonparallel planes, and convex ring-torus sectors closed by two meridian
disks. These parts previously returned `field.import-domain`. Face names,
membership and analytic volume are available through the existing interfaces.
The [M21 note](../M21-CLIPPED-IMPORT.md) gives the exact admission rules and the
remaining refusals.

Existing primitives and one-axis revolutions keep their earlier path. The new
reader runs only after a revolution-domain refusal. It checks complete shared
circle boundaries, orientation and incidence before retaining the intersection
of analytic supports. The field does not retain imported topology.

Unsupported spherical cuts report the separation gap and squared margin.
Unsupported torus sectors report their cross-plane signs. General spline faces
keep `field.import-spline` behind M10. M16 solver queries on imported parts are
parked and are not extended here.

No wire format, tolerance, precision policy or work allowance changes. Raw
evidence retained privately.
