# F7 solver query adapter

Adds the opt-in `catallax_field::solver` interface described in
[Solver queries](../SOLVER-QUERIES.md). `Definition` owns a revision-bound native
plan and explicit material assignment; `Queries` exposes point/cell membership,
finite-cell integration, true-distance and named-face services, and material
membership. The grid uses transferred operation allowances and fixed ordered
jobs, with all available cores by default.

No existing public method changes. The unchanged private mass axis-box detector
becomes crate-visible for exact box integration; mass decisions and formulas are
unchanged. No print-cell selection, saved answers, allowances, tolerance,
constitutive solver, proof consumer or export decision changes. Existing FEA and
reconstruction examples remain untouched. Named patch integration and fast exact
closest construction support a single native axis box; full M16/M20 coverage is
not claimed.
