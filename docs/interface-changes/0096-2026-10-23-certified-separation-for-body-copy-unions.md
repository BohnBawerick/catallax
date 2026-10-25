## 2026-10-23: certified separation for body-copy unions

Pieces: features, measure
Change: feature evaluation reuses the existing certified `bounding_box` query for a pattern's seed and for mirrors that retain the original. Exact affine enclosure propagation and strict box separation can prove copies disjoint, in which case union retains their immutable lumps directly. Other copies still use the existing Boolean engine.
Why: disjoint trimmed bodies need no intersection of their unbounded carriers. Avoiding those intersections removes needless branch-tracing refusals and request-work exhaustion.
Migration: no query signature, geometric tolerance, work allowance or existing Boolean operation changes. Enclosures prove separation only; loose or touching boxes never prove a geometric contact or move a face.

