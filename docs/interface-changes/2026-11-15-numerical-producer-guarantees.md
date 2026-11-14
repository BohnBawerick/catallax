# Numerical producer guarantees

`build_numerical` and numerical continuation return `NumericalBody` with
`NumericalEvidence`. Their public summaries report `exact: false`, including
when every sampled residual is zero. They do not produce `BoundedEvidence` or
grant exact checkpoint or STEP export authority.

This is a guarantee about those producers. The existing generic `Answer` and
`BoundedEvidence` constructors trust their callers. An explicit caller can
mislabel a numerical body as exact or a sampled value as a complete bound.
Preventing that misuse requires a separate producer-owned result-authority
change across the shared contracts and their consumers. N01 retains those
contracts. Its compile-fail examples check the absence of implicit conversions,
not all possible explicit uses of the generic wrappers.
