## 2026-11-12: Exact STEP spline admission through existing certificates

Pieces: 2, 3, 8

Change: STEP solid construction tries spline support, finite endpoint, embedding, incidence and boundary certificates before the existing self-intersection and topology checks. Unsupported spline geometry names the failed certification step in the existing `Unsupported` payload.

Why: the blanket spline gate rejected definitions before the kernel could attempt its available proofs.

Migration: callers keep the existing signatures and serialized `ImportRefusal` variants. Treat an unsupported spline as a failed named proof step, rather than assuming that every spline solid is unsupported. A returned body still requires every solid check to pass. Source previews remain unproved and are not certificate inputs.

This extends the construction limit recorded in the [source-preview fragment](2026-11-07-step-spline-source-preview.md). Source recognition, exact decimal definitions, canonical source hashing and the original file bytes retain their existing contracts.

A spline surface must pass `RationalPatch`'s whole-domain embedding proof. Its physical face must then pass `RationalPatchFace`'s complete four-row identity, endpoint and material-relative orientation checks. Arbitrary trims, periodic supports and patches without a global embedding certificate retain named refusals. Declared STEP closure and self-intersection flags do not supply proof.

A spline edge needs unique finite endpoint preimages, consistent parameter direction and a strictly monotone projection over every knot span. Collinear certified arcs can use an exact line carrier inside the checked body. Other arcs keep their exact spline definition. Incidence on planes, cylinders, spheres and ring tori is a homogeneous polynomial identity over every span crossed by the finite edge. Cone incidence retains a refusal because the implicit equation alone does not prove nappe membership.

Spline faces use the complete positive-weight surface control hull for broad-phase exclusion, including interior controls. Certified boundary facts can serve the existing M10 intersection engine, but a missing downstream intersection or trim query still returns `SelfIntersectionNotVerified`. A successful support or boundary proof alone never admits a solid.

Spline admission performs no fit or repair and adds no nonzero local tolerance. The 1 nm modelling tolerance, approximation caps, M24 repair rules and caller work allowances retain their existing values. Exhausted operation work returns `ResourceLimit`.
