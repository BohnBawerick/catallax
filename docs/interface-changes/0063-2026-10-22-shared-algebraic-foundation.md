## 2026-10-22: shared algebraic foundation

Pieces: 1 to 7, field, display mesh, two-kernel comparison and cli consumers; judge formats and kind tags unchanged
Change: implement the algebraic-definition decision across the existing consumers. `RealPoly` stores ordered `Real` coefficients, with certified degree, equality, resultants and univariate root selection; rational polynomials retain the existing `MPoly` backend. `sin_cos_degrees` constructs certified algebraic sine/cosine values from exact rational degrees. `Decider::rational_value` optionally proves an exact rational value without replacing its expression definition; `None` is inconclusive, not a proof of irrationality or definedness.
Why: frames, poses, analytic dimensions and their support equations must describe the same exact geometry in both kernels, including irrational rigid motions and oblique circular extrusion.
Migration: use `Point3`, `Vector3` and `Real` for stored definitions, including `PointDef::Explicit` and `CurveForm` directions. Analytic constructors and `Pose::new` take a final `&mut Decider`; `Rotation::from_real_quaternion` adds checked algebraic rotation input. The rational quaternion constructor remains available. `matrix()` exposes `Real` entries; `rational_matrix()` and `rational_components()` borrow exact rational entries when available. Definitions containing `Real` have no mathematical equality/hash traits; use certified geometry comparisons. Construction identity can prove cache reuse only.

`Pose::apply_rational_point` and `apply_rational_vector` preserve native rational-input arithmetic without temporary expression leaves. They return `None` when the stored frame does not expose rational entries; they do not round or decide whether an algebraic expression happens to be rational. The existing recipe feature constructors use this path; runtime field placements and ray samples remain algebraic.

`Pose::from_rational` accepts a validated rotation and a native rational translation without a numeric proof. Rational recipe sampling uses this constructor and the native point transform.

Charts take the caller's decision context and return an outer `Result`, retaining their former inner optional result where applicable. Implicits and curve equations return `RealPoly`; consumers use `rational()` or `into_rational()` for the borrowed or consuming rational fast path, and certified operations otherwise. No rational conversion rounds an algebraic definition.

`EllipticCylinder` adds its normalized transverse implicit, periodic/linear chart, plane intersections, measurement over the actual oblique trims, and STEP `SURFACE_OF_LINEAR_EXTRUSION` over an `ELLIPSE`. Equal semi-axes reduce to a circular cylinder before construction. New-kind display meshes and new recipe coverage remain separate work.

Field plans retain algebraic support polynomials, gradients, ray germs and rigid placements. `compile_checked`, `proven_checked`, `face_query_checked` and `transformed` accept the shared decider. Convenience queries use the default policy and return an explicitly indeterminate, unbounded result when a numeric proof fails; checked queries preserve the numeric refusal. Rigid transforms preserve the existing metric formulas; no affine-stretch distance shortcut is introduced.

The comparison consumer uses `RealPoly` support equations and gradients, certifies the first nonzero normalisation coefficient, and compares supports through the shared decider. Rational equations retain the original `MPoly` arithmetic. Rational-only witness and feature-candidate APIs do not round algebraic placements: an unavailable native candidate frame returns `field.sample-algebraic-placement`, a sampling failure rather than a geometry refusal or an agreement claim.

Algebraic field coefficient, metric and ray-germ enclosures use `Decider::enclose_budgeted_at`, preserving the caller's precision and work limits. The existing `enclose_at` remains the explicitly unbudgeted, known-defined working helper. Native rational enclosures keep their original arithmetic path.

Intersection-plan caches retain input definition nodes. Their keys therefore use `Surface::shares_definition_nodes`; value-equal independent definitions cannot inherit another query's expression-DAG aliases and change charged proof work. `Real::same_node` is a storage-identity query, not a mathematical equality test.

This supersedes the rational-only primitive-storage statement in the preceding conic-dimensions entry. Question/manifest version 2, supported CLI kind tags, numeric limits and deadlines are unchanged. Elliptic cylinders and open conics are not exposed through new or remapped CLI tags.

