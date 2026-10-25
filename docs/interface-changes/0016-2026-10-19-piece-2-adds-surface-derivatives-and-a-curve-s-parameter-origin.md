## 2026-10-19: piece 2 adds surface derivatives and a curve's parameter origin

Pieces: 2, 4, 6, 7
Change: `Surface::derivatives(&self, u: &Ball, v: &Ball) -> [[Ball; 3]; 2]` encloses `[dS/du, dS/dv]`, and `Curve::parameter_origin(&self, store: &GeomStore) -> Point3` gives the exact point at `t = 0` (a circle's or ellipse's point in its reference direction, an intersection curve's witness).
Why: piece 6 integrates over surface parametrisations and needs their derivatives; pieces 4 and 7 start closed curves and place ring-edge vertices at the parameter origin, which needs one exact definition.
Migration: piece 2 creator implements both. Pieces 4, 6 and 7 use them instead of deriving their own.

