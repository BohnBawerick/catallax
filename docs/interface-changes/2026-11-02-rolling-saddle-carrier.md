## 2026-11-02: fitted ball carrier between two nonparallel cylinders

Pieces: shapes
Change: `catallax_shapes` gains `RollingSaddleDefinition`, `RollingSaddle`, `BallSide` and `ROLLING_SADDLE_TARGET_NM`. `RollingSaddleDefinition::new` takes a wide and a narrow cylinder, a ball radius, the ball's side of each and the ring to follow, and `section` gives the exact ball centre and both contacts at a rational chart parameter. `RollingSaddle::fit` builds a closed rational carrier of degree (9, 2) with rails exactly on both cylinders, a recomputed two-sided distance bound and a regularity proof, and `rails`, `bound`, `span_count` and `transformed` read or move it. Cases outside the regular ring refuse with `shapes.rolling-saddle`, and a span the proof cannot resolve refuses with `numbers.geometric-gap-unproven`. This is not yet a `SplineSupport` variant, so no recipe reaches it.
Why: this is the first piece of the nonparallel cylinder slice. The centre path of such a ball is an elliptic curve, so no exact rational surface exists, and this is the kernel's first fit for a closed ring between two curved faces.
Migration: none. Nothing that compiled before changes.
