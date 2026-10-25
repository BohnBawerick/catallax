## 2026-10-26: exact annular torus first-hit limits

Pieces: features, field

Change: `extrude_to` and `extrude_next` accept coaxial annular profiles against a stock torus when one entry or exit branch is certified over the complete profile. Both kernels use exact circular sections and finite construction columns. Authored direction scale and sign, rigid placement, profile holes and extrusion-owned names are retained. Oblique or offset torus sections and spline intersection pairs remain named coverage gaps. The quadratic cap reduction now refuses higher-degree ray terms instead of omitting them.

Why: planar limits and rectangular stock-quadric caps did not cover annular features stopping at a toroidal wall. The kernel already constructs the required coaxial cylinder/torus and normal-plane/torus circles.

Migration: none. Recipe versions, work allowances, modelling tolerance and existing planar and stock-quadric definitions are unchanged. Public test parts cover each curved surface family; represented spline limits explicitly retain the unsupported-pair refusal.

