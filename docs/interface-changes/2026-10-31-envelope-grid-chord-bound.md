## 2026-10-31: envelope display grid follows the carrier's curvature

Pieces: shapes, mesh
Change: `RollingEnvelope::display_grid` now bounds a cell by the Bernstein coefficients of `P - B W` over the smallest weight instead of comparing rational controls with their Greville points. The P5 radius 2 grid at 1/200 mm has 1024 cells instead of 4096, and its display mesh at the default 1/100 mm has 40,776 triangles instead of 57,096.
Why: the old comparison added a term from the weights' own slope, so it refined past the carrier's real curvature. Every application request on a fitted body pays for this grid, and SolidDraft's screen view spent minutes on it.
Migration: none. Stored bounds stay two-way and at or below the requested tolerance. A consumer that read a grid or display mesh size sees smaller ones. Certified exchange, STEP and volume do not read this grid.
