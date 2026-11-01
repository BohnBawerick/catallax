## 2026-11-02: exact straight-path sweeps

Pieces: features, application, comparison, measurement, display, STEP

Change: The B-rep recipe evaluator accepts an empty path or one nonzero rational line segment, and evaluates `sweep` along it. Both `fixed` and `rotation_minimising` have the same identity rotation on one straight segment. Recipe versions 3 and 4 keep their existing schema. No public Rust signature changes.

Why: This is the smallest analytic sweep with complete existing predicate, integration, display and exchange consumers. It does not require a general spline carrier or an approximate path.

For start `a`, end `b` and profile region `R`, the motion is exactly `q + t (b-a)`, with `q` in `R` and `t` in `[0,1]`. The profile stays at its authored placement, not at `a`. A nonzero normal component makes this affine map injective over the complete profile domain, so the existing extrusion shell is exactly the swept solid, including nonconvex boundaries and holes. Its model displacement is zero. No square root or unit-direction approximation changes the displacement. A profile-plane motion has zero three-dimensional volume and regularises to the empty body, as does an empty path or profile.

The immutable recipe retains the path and sweep definitions. Incremental checkpoints retain the checked path value. Caps use the sweep step's `start` and `end`; side origins include both the profile `segment` and the authored path segment `path`. Shared physical edges keep their names. Existing planar, circular-cylinder and elliptic-cylinder consumers operate on the exact analytic reductions. Display and STEP each retain their separate requested bounds; neither bound replaces the zero model-displacement proof.

A zero-length line refuses at the path step with `path.zero-length`. Curved paths and paths with more than one segment refuse there with `features.unsupported-path`, a coverage gap, not impossible geometry. Lofts, spline and helix path syntax, turns and general derived carriers remain unimplemented. The field compiler does not implement recipe paths or sweeps; comparison reports that missing field coverage rather than claiming independence.

Migration: No existing extrusion caller changes. Applications may now emit the existing `path` and `sweep` recipe steps, preserve the authored profile placement and use the endpoint difference as the intended motion. Do not normalise that difference into a unit vector. Use the ordinary `catallax_app::build` or incremental build, picking names, requested measurement and `catallax_step::write_step` paths. The `soliddraft-recipe/1` converter has no new operation; SolidDraft must emit Catallax recipe steps directly or add its own conversion before exposing the feature. Only enable the one-line path domain, and keep unsupported paths visibly refused.
