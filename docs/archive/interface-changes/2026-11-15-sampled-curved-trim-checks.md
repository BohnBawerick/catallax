Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Sampled curved trim checks

Numerical admission checks planar spline loops through an eight-chord walk per
edge. The resulting polygon can be simple and contain its hole while the actual
curve leaves the outer boundary between sample nodes.

For these projected spline walks, `trim.sampled-hole-contained-disjoint` now
reports the polygon check. The producer emits no unqualified
`trim.hole-contained-disjoint` entry. Continuous containment and disjointness
remain unresolved. `NumericalCheck` records only checks actually performed;
matching nodes do not establish the intervening geometry.

Linear parameter trims retain their existing unqualified polygon check. The
discrete boundary-endpoint check and pointwise distance samples keep their scope.
Tolerance, sampling, geometry admission and result category are unchanged.

The regression uses a native spline hole edge that matches all nine sampled
locations but reaches `x = -1` between nodes, beyond its outer rectangle. Admission
still reports its actual numerical checks and retains no continuous hole claim.
