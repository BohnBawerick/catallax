## 2026-10-21: regular higher-order radial contacts

Piece: 5
Change: after tangent and curvature ties, regular germs may use exact implicit Taylor coefficients through order eight. Both use the same transverse linear parameter, so unequal tangent speeds cannot change the order. Singular or unresolved germs and shared components still refuse.
Why: the remaining contact case has a common sphere and internally tangent cylinders of radii three and two. With transverse parameter `z`, both branches have `y = z^2/4 + z^4/64 + ...`, but their `x-3` coefficients are `-z^4/96` and `-z^4/64`. Second-order predicates cannot distinguish them.
Migration: no public Rust signature or default budget changes. The internal radial sorter receives the model and each germ's source. Tests derive third- and fourth-order conic coefficients independently, vary tangent speed, reverse comparison order and retain refusal for duplicate germs.

