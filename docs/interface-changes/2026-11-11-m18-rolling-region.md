# Exact finite lower rolling-ball region (partial M18)

Adds `catallax_field::rolling_region::{Region,Station,Relation,Record,Batch}`.
This is an opt-in explicit mathematical region with strict membership/refusal,
exact contacts and an ordered parallel query service. It does not compile a
recipe fillet or establish P5 acceptance. See [the set, proof obligations and
limits](../M18-ROLLING-REGION.md).

Adds read-only `Decider::deadline()` so transferred parallel numeric jobs inherit
the caller's absolute deadline as well as its existing policy. No defaults or
allowances change. `Region::classify_batch` consumes one caller-owned allowance,
partitions by input index independently of thread count, and returns complete
ordered per-point results, work and exhaustion. Default execution uses all
available cores up to input count; recipe compilation remains serial.

No named-face, distance, mass, print, export or existing proof route delegates
to this new service. Existing recipe refusals are retained.
