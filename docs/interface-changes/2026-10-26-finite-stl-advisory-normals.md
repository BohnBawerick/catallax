## 2026-10-26: finite STL advisory normals

Pieces: field print writer and B-rep print adapter.

Change: `catallax_field::printing::Mesh::write_stl` preserves every previously finite advisory-normal byte. If binary64 subtraction and cross-product cancellation would produce a non-finite normal, it computes the cross product from the exact encoded binary32 coordinates before normalization. Non-finite coordinates and genuinely degenerate facets return `std::io::ErrorKind::InvalidData` instead of emitting NaN normal values.

Why: a valid finite binary32 triangle can lose its nonzero normal during binary64 edge subtraction. The triangle with vertices `(2^60, 2^60, 0)`, `(1, 0, 0)` and `(0, 0, 0)` has exact normal `(0, 0, -2^60)`, although the old floating-point calculation returned zero. It is a face of a closed tetrahedron with fourth vertex `(0, 0, 1)`.

Migration: no signature, request, tolerance, allowance, metadata or finite field-print byte changes. B-rep print qualification remains authoritative for facet winding, topology and material orientation. The shared writer now also guarantees finite advisory-normal coordinates. The CLI removes its print artifacts on writer failure. This supplements the [B-rep print contract](2026-10-26-certified-brep-print.md).
