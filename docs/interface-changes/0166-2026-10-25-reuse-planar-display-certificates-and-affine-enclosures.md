## 2026-10-25: reuse planar display certificates and affine enclosures

Pieces: display mesh, application benchmark

Change: planar facet quality checks use a cached interval unit normal, with the existing multiprecision certificate retained when the filter is inconclusive. Shared-rail ear centroids enclose the same exact affine construction from cached endpoint balls instead of expanding exact trace coordinates. Requested edge-length refinement still retains the exact centroid it consumes.

Why: fine cap meshes repeated multiprecision angle constants per triangle, and long rational trim coordinates made repeated exact centroid construction expensive.

Migration: none. Public mesh options, user tolerances, resource limits, shared edge polylines, whole-facet bounds and exact geometry are unchanged. Rounded interior coordinates can change within their certified bound. The display benchmark accepts optional chord and angle arguments and records per-face bounds; its default invocation is unchanged.

