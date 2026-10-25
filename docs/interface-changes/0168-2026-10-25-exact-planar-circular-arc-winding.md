## 2026-10-25: exact planar circular-arc winding

Pieces: 5, application

Change: planar loops mixing lines and circular arcs use signed chord crossings
plus the exact circular-segment winding correction. The correction uses only
disk and chord-half-plane signs; it constructs no new intersection geometry.
Minor and major arcs, reversed loops, holes and seamed full circles retain their
closed-trim semantics. Selected-arc and finite-chord queries keep the existing
boundary fallback, and unsupported carriers keep the existing probe/trace path.

Why: rounded planar caps missed the straight-polygon and complete-ring paths.
Repeated membership queries constructed auxiliary surface probes and spent the
application's aggregate arithmetic allowance on those intersections.

Migration: none. Public APIs, ordered section carriers and points, naming,
work allowances, geometric tolerances and display options are unchanged.
The public rounded base/flange recipe now completes the default full application
build; its regression checks the named inner joint edge at exactly `(0, -35, 10)`.
Timing and comparison tooling remains outside part construction.

Prior art: [Franklin's crossing convention](https://wrfranklin.org/Research/Short_Notes/pnpoly.html)
keeps ray endpoints half-open. The [circular segment](https://mathworld.wolfram.com/CircularSegment.html)
bounded by an arc and chord gives the algebraic winding correction. No source
code is copied and no angular or sampled approximation is introduced.

