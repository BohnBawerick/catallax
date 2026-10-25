## 2026-10-21: retained exact probe incidence

Piece: 5
Change: probe-plane constructions record the incidence they establish, and rational face probes retain the exact side-query proof. For an unmarked surface pair with one curve, existing incidence on both surfaces proves query incidence on that curve. Pairs with multiple components retain the branch-membership predicate. Construction provenance separates auxiliary probe geometry from intrinsic seed supports, preventing later probes from recursively following earlier probe sections while retaining the facts themselves.
Why: rebuilding an already-known conic in a different algebraic coordinate representation forced another zero elimination and exhausted the decision. The missing information was the existing incidence proof, not a geometric tolerance.
Migration: no public API or budget changes. Derived facts retain exact provenance. The two-generator meridian regression compares public curve membership with independent line-distance predicates, so support membership cannot be mistaken for membership of the wrong component.

