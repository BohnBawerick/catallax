## 2026-10-26: bounded circular chamfer strips and shared revolution compatibility

Pieces: features, field, recipe consumers, STEP regression tests.

Change: Convex finite circular plane/cylinder chamfers construct exact conical strips, contact arcs and shared planar end boundaries. Finite concave circular chamfers return `chamfer.unsupported-concave-arc` in both kernels. Existing whole-ring and straight chamfer domains retain their contracts. Field consumers add bounded exact straight-loop revolution and validated radial-sector membership paths; unsupported profiles keep their existing fallback.

Why: Endpoint trimming previously required straight rails even for circular edges whose exact conical support was already available. This corrects the trimmed-circular implementation gap recorded in [the elementary chamfer entry](0135-2026-10-23-exact-elementary-chamfer-construction.md); general overflow and corner continuation remain outside this slice.

Migration: No recipe schema, public signatures, work allowances, modelling tolerance or output accuracy changes. Consumers must check the full requested result. This draft remains parked: the original convex native display exceeds the unchanged 1M scoreboard allowance, and the chord fixture exhausts the separate 100M verifier allowance. Finite concave support needs separate end-domain and native-consumer qualification. STEP read-back and sampled agreement do not certify whole-result deviation.

Shared-hook decision: PR93 landed first. This second-lander integration retains its `Box<Angle>` representation and exact Angle API while preserving partial-sector eligibility in `convex_revolve`. The full-turn-only eligibility guard is removed; angular classification uses `Angle::is` and `Angle::location`. Regressions exercise the actual public chamfer section at a partial turn and full-turn recognition with irrational axis normalization. See [the revolution contract](2026-10-25-algebraic-revolution-frames-and-radii.md).

Membership reconciliation: Preserve P3's corner/ring classifiers, elementary fillet inequalities and interval exclusion proofs. The chamfer sector cache is confined to the existing extrusion path and applies only when the existing corner classifier does not apply. Both lanes' independent membership, face and ray-germ regressions remain enabled. See [P3 native display](2026-10-25-p3-native-field-display.md). This integration does not reopen chamfer optimization or clear readiness, replay or qualification gates.
