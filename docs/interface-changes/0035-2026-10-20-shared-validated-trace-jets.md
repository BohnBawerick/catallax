## 2026-10-20: shared validated trace jets

Pieces: 4, 6
Change: trace clones share a bounded recent cache of Taylor enclosures keyed by exact parameter midpoint, radius, working precision and order. Constant enclosures are cached separately for each requested precision. Measurement traces each edge once per scalar-property call, and shares the traces across volume, area and centroid within `mass_properties`. Each property retains its fresh arithmetic budget.
Why: full-question profiles placed most measurement time in repeated trace Taylor evaluation. Reusing an already certified enclosure removes that computation without changing the enclosure or its proof. Precision-specific constant caching also prevents call order from changing an enclosure.
Migration: no public signature, requested accuracy, refusal code or budget default changes. Cached requests still check the deadline. Cache eviction changes cost only.

