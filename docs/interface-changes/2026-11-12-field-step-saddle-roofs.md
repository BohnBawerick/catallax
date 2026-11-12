## 2026-11-12: the field builds imported saddle roofs, with the quadric of a curved bilinear patch

Pieces: field, compare

Change: `Plan::compile_with_files` builds an imported STEP part that is a saddle roof: a curved bilinear equal-weight spline patch, a floor plane and four wall planes that contain the patch's twist direction. The field holds it as the intersection of the six supports, the patch by its exact quadric, under the file's face names, and `field::mass::volume` gives its exact volume. `catallax_field::support_polynomial` returns the support equation of a bilinear spline patch (its plane when flat, its quadric when curved), `patch_witness` a point of a curved patch near a given point, and the field's candidate sampling covers the patch. `catallax_compare` uses them to check a curved spline face of a mesh against the field.

Why: M21 asks for imported parts with spline faces to build in the field, with the two kernels agreeing on request. This extends [flat spline patches](2026-11-12-field-step-flat-spline-patches.md) to the first curved patch the shared importer admits.

Details:

- A bilinear patch is `a + u b + v c + u v d`. In the basis `b`, `c`, `d` with coordinates `alpha`, `beta`, `gamma` it is `gamma = alpha beta`, so its support is the quadric `T^2 (gamma - alpha beta)` with `T` the triple product of `b`, `c`, `d`, signed along the patch normal. The planar walls are `alpha = 0, 1` and `beta = 0, 1`, and the floor is `gamma = p + q alpha + r beta`.
- The part is accepted only after an exact check: six faces, one curved patch, walls on exactly those four planes with outward normals, a floor strictly on one side of the patch at all four corners (a bilinear difference takes its extremes at the corners), every face meeting at the eight corners these supports give with straight edges, and a point between roof and floor strictly inside every halfspace. Then the imported solid and the intersection are bounded by the same faces on the same supports and are the same solid. Anything else, including a curved patch in a solid of any other shape, is `field.import-spline` with the reason.
- A reversed face sense turns the field's patch over (its parameters swap), so the support's gradient always points out of the material.
- The volume is `|T|` times the mean gap between patch and floor over the unit square, a rational number.
- The quadric's distance residual is a lower bound, not a distance: `2|F| / (G + sqrt(G^2 + 2 K |F|))` with `G` the gradient length and `K` a bound of the Hessian norm, derived from `|F(p)| <= G d + K d^2 / 2`. It never exceeds the distance to the whole quadric.
- Higher-degree patches and rational patches (unequal weights) stay refused until the field side of M10.

No existing answer changes: every part this admits was refused before. `field::mass::volume` keeps its results for every other leaf. No tolerance, saved replay answer or work allowance changes. Raw evidence retained privately.

Migration: Callers that treated `field.import-spline` as final for a curved patch now receive a built field for a saddle roof. Code that matches exhaustively on the field's internal support kinds is unaffected; the new kind is private.
