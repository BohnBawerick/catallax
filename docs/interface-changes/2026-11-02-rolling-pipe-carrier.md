## 2026-11-02: exact rolling pipe between a plane and an oblique cylinder

Pieces: shapes
Change: `catallax_shapes::RollingPipe` is a new public type. It holds the exact ball-centre and contact ellipses of a closed regular contact between a plane and a cylinder that meet obliquely, and an exact rational degree (4, 2) `BSplineSurface` carrier over `[0, 4] x [0, 1]`. Row `v = 0` is the plane contact and row `v = 1` is the wall contact. `RollingPipe::new` refuses with `shapes.rolling-pipe` when the plane is not strictly oblique, when the cylinder offset reaches its axis, and when the ball could fold at the tightest curve of its centre path.
Why: general fillets need a first carrier that is not elementary. The centre path here is an ellipse and its normal plane is spanned by two rational unit rays, so the surface has a closed rational form and no fit or fitted error.
Migration: none. Nothing calls this type yet. `FilletDefinition` and the recipe `fillet` operation still refuse these edges, and `SplineSupport` has no new variant, so no consumer needs a new match arm.
