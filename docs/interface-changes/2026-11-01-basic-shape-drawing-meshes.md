## 2026-11-01: full periodic faces and tilted cylinder end rings

Pieces: mesh, feature integration

Change: loopless spheres and tori that exhaust the general periodic triangulator can use a structured periodic lattice when no edge-length limit is requested. Seam indices wrap without temporary full-period copies, and a sphere's pole rows collapse to fans. Floating estimates choose sample counts only. The mesh retains the existing rounded-position residual, Hessian chord bound and per-facet normal certificate. Counts and refinement remain bounded by the caller's unchanged vertex cap. The general route still runs first, so its successful outputs do not change.

The general triangulator also rejects an oversized untrimmed periodic seed grid before allocating its chart copies. This changes when the existing vertex refusal is detected, not which grid can finish.

A cylinder with two analytic rings, including an ellipse, can use the existing certified angular merge stitch. This is the same stitch used by fitted cylinder rails. Its winding, orientation, chord, normal and optional edge-length checks remain in force. Cylinders with only circular rails retain their previous route. Successful meshes with analytic elliptical rails can change; their new stitch reports the achieved certificates rather than retaining a triangulation grown by centroid refinement.

Why: the drawing chord and 2-degree target exhausted the grid or temporary-copy budget on an R1 full-revolve sphere and R3/r1 full-revolve torus. An R1 cylinder cut through z=5 at 28 or 53 degrees seeded one axial cell and refined toward a million chart vertices instead of stitching its two rings.

Migration: no public API, tolerance, vertex cap, work allowance, construction geometry, STEP output or application setting changes. `crates/catallax-features/tests/drawing_basic_shapes.rs` covers the supplied full-revolve profiles and the original height-10 rod cutter, plus the 53-degree variant, at native display, automatic viewport and drawing controls. It checks closed opposite edge uses, enclosed material volume, vertex limits and returned chord and angular certificates. Full-suite and replay evidence is produced by CI.
