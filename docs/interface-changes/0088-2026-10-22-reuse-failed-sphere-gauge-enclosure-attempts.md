## 2026-10-22: reuse failed sphere-gauge enclosure attempts

Piece: 6
Change: grouped volume, area and centroid measurement retain request-local sphere-gauge domain failures by face, subdivision depth, rectangle/Taylor mode, precision and gauge. A later property skips an attempt only when all six gauges already failed the identical area-gauge denominator check and every boundary arc uses the same trace integration path. Property-specific analytic arcs and second moments do not use this cache.
Why: centroid measurement repeated the same failed domain enclosures that volume and area had already resolved by subdivision. Those failed attempts consumed its measurement allowance before the first finite centroid enclosure.
Migration: no public signature, refusal code, requested accuracy, default allowance or deadline changes. Successful panel arithmetic and retry ordering are unchanged. Individual failed gauges are not skipped inside a successful attempt, so its quadrature-dominance flag remains unchanged. A failed enclosure is never treated as proof that a geometric point or region is absent.

