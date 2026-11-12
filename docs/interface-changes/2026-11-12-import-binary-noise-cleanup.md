## 2026-11-12: recorded binary-noise cleanup for curved imports

Pieces: 8, features (import steps), STEP and IGES consumers

Change: when the exact reading of a STEP or IGES file refuses for incidence
(`not_solid` with "exactly"), a rotation that is not exactly orthonormal, or a
self-intersection, `apply_source` tries a cleaned reading of the same bytes. Each
decimal with at least ten significant digits becomes the simplest rational within
a relative `1e-12`, when its denominator is at most 10,000. The cleaned file must
pass every exact check of an ordinary import. Its whole-file movement is proven
(see [piece 8](../pieces/08-step-import-repair.md#binary-noise-cleanup)) and must
fit the 0.01 mm automatic cap and a tenth of the file's conservative feature
size. Every face is then recorded in `face_tolerances` with that bound, and the
cleaned file is the derivative. `catallax_import` gains `rebuild_source`;
`rebuild` is its STEP form. Files with spline data are not cleaned.

Why: PLAN M24 asks for curved repairs with recorded local tolerances up to
10 µm. The R8 census's `not_solid` refusals include vertices and edges that miss
their curves and faces by binary rounding. OpenCASCADE's cylinder, cone, sphere
and torus IGES files all refused for this reason.

Coverage, raw evidence retained privately: OpenCASCADE 8.0 IGES cylinder, cone,
sphere and torus now import through the cleanup with their analytic volumes, a
proven per-face bound below 1e-9 mm, and a derivative that `rebuild_source`
reproduces. A STEP cylinder whose ring vertices carry `-7.3e-16` imports the
same way, while a deliberate 1e-6 offset written with eleven digits and a short
`1E-15` are not snapped and keep their `not_solid` refusal. The exact reading of
every one of these files still refuses; only the repair entry cleans. The R8
census was not rerun here.

Migration: none for callers of `apply`, `rebuild` or the recipe steps; they may
now return bodies, records and derivatives for files they refused before. Files
that imported before are read exactly as before. No tolerance, cap or work
allowance changes.

Known gaps: noise around values that are not simple rationals, such as rounded
`sqrt(2)/2` frames, still refuses; that needs vertex snapping onto the exact
curves, not decimal cleanup.
