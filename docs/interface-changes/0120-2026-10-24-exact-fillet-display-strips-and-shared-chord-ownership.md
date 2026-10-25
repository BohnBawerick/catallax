## 2026-10-24: exact fillet display strips and shared-chord ownership

Pieces: display mesh, application regression coverage

Change: finite cylindrical fillets use the same certified ruled-strip triangulation as elliptic cylinders, retaining every shared end-curve sample. Toroidal bands lift the minor angle into a nonperiodic interval selected by their oriented loops. Certified frame and angular evaluations are cached locally for that band; grid dimensions follow the weighted curvature terms rather than the cone's aspect heuristic. The existing half-tolerance grid reserve and mandatory triangle/edge certificates are unchanged.

Curved-face triangulation also records interior chords between boundary vertices in face order. A later face refines an unconstrained chord already owned by another face. Shared B-rep polylines remain fixed. This prevents adjacent cylinder fillets from publishing opposite copies of a boundary-only triangle and giving their interior diagonal four uses.

Why: P3 previously refused mesh closure at adjacent cylinder strips. Exact convex torus fillets displayed but traversed unnecessary periodic copies and repeated certified angular evaluations. The native application regression now exercises P3's four cylinder blends and torus blend, face picking roles, shared-edge incidence and radius-edit names.

Migration: no exact model, fillet construction, API, naming, tolerance, work allowance, manifest or STEP-option changes. General curved-boundary acceleration remains separate work.

