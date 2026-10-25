## 2026-10-22: proven centroidal inertia on request

Pieces: measurement, application and question measurement publication
Change: `MassProperties` includes `inertia: Option<InertiaTensor>`, a symmetric 3-by-3 enclosure matrix about the centroid at unit density. Empty bodies have no centroid or tensor. `inertia` requests only the tensor; `centroid` requests only the centre. `volume_area_centroid` preserves the original grouped measurement without second-moment work. `ApplicationBody::mass_properties` and `ApplicationBody::inertia` borrow the built result, so a refusal does not discard it.
Why: applications need proven inertia for mass-property reporting and equivalent B-rep bodies, rather than mesh estimates or guessed values.
Migration: callers of `mass_properties` receive the additional measurement and at most one additional width mark. Call `volume_area_centroid` to retain the original property set and budgets. Questions opt in with `measure.properties: ["inertia"]`; the matrix uses model axes and length units to the fifth power, with negative products of inertia off the diagonal. Default question measurements, other verification answer formats and STEP definitions are unchanged.

