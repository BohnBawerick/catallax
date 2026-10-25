## 2026-10-19: piece 2 fixes its exact forms and validation order

Pieces: 2, 3, 4
Change: `Surface::implicit` and `Curve::equations` return exactly the polynomials of the piece 2 tables, unscaled and in the stated order; a line's equations are the non-zero components of `(x - origin) x dir` (two or three) instead of "two planes". An ellipse now needs `0 < minor_sq < major_sq`; equal semi-axes are refused with `shapes.ellipse-axes` and must be a `Circle`. Constructors check zero vectors, then scalar ranges, then perpendicularity, and return the first failure. `parameters_of` uses the stated formulas, with `u` in `[0, 2 pi)` so its jump is on the seam. `Rotation` equality compares the matrix.
Why: blind tests need one exact answer per call, and one representation per circle keeps `same_surface` and piece 4's closed forms unambiguous.
Migration: piece 4 builds an ellipse only when its axes differ and a circle otherwise. Piece 3 builders validate before calling piece 2 constructors, so they meet none of these refusals.

