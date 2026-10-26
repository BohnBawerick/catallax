## 2026-10-23: spline carriers integrate with trimmed analytic distances and drill cells

Pieces: shapes, Boolean filters, features and field
Change: analytic `distance::Domain::new` and `distance::project` return `shapes.distance-uncovered` for finite B-spline supports rather than assuming a global analytic frame. The bounded-face exclusion filter uses an implicit equation only when one exists. Definition-native drill and capped-sweep cells propagate fallible constraint construction through their exact isometry API.
Why: the independently landed analytic-distance and feature consumers must retain their coverage guarantees when the genuine spline kind is present.
Migration: generic trimmed analytic distance charts do not yet cover represented spline faces. Use the separately certified original-target prism field for its nominal true-distance contract; do not substitute that definition for represented-face distance. An unavailable analytic exclusion proof falls through to the existing intersection path and never proves disjointness. Existing analytic definitions, numeric allowances and output accuracies are unchanged.

