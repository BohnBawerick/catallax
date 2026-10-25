## 2026-10-24: analytic columns for certified field moments

Pieces: field, CLI stand-in exchange
Change: Constant axis-aligned cubic grids with Boolean axis boxes and z-parallel cylinders can use certified analytic-column integration. `mass::Options::max_cells` and `Properties::cells` count classified boxes or analytic intervals, including initial partitions. The rational estimate is a moment proposal, with exact surrogate and exchange enclosure checks still required.
Why: Voxel boundary uncertainty alone does not reach strict inertia widths efficiently. Inner antiderivatives and interval Taylor remainders reduce the integration dimension without using print geometry.
Migration: No caller or format change. Requested Accuracy values, defaults, work allowances and the strict 1e-9 stand-in inertia comparison are unchanged. Other field definitions retain box integration. Measurements remain explicit queries, not part-build work.


