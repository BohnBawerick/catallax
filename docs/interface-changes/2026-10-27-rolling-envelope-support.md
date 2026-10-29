## 2026-10-27: checked rolling-ball envelope spline support

Pieces: 2 (catallax-shapes), and every piece that matches `SplineSupport`
Change: `SplineSupport::RollingEnvelope` adds a checked finite rational B-spline family for a constant-radius ball that rests on a cylinder and on a parallel circular mouth edge. It comes with `RollingEnvelopeDefinition`, `RollingEnvelope`, `ROLLING_ENVELOPE_TARGET_NM` and the refusal code `shapes.rolling-envelope`. The new `SplineSupport::param_kinds` reports `[Finite, Finite]` for this family, and `Surface::domain` now delegates to it.
Why: the P5 housing junction ball rests on each countersink mouth between its plane-to-edge transitions. That envelope is not rational, so it is represented by a fitted carrier:
- The cylinder rail, mouth arc and transition meridians are built exactly on their supports.
- The carrier's two-sided position bound against the exact envelope is recomputed from interval Taylor jets, targeting 100 nm, with a 1 um cap.
Migration: matches on the non-exhaustive `SplineSupport` keep their wildcard arms. A consumer that needs this family adds its own arm. No existing elliptic-offset behaviour changes.
