# Torus interior probe order

STEP chart preparation interleaves the two existing rational ray families
when probing through a partial torus face. It retains the four axis rays,
every later candidate and all crossing, orientation and boundary checks.
Chart cuts keep their original candidate order and topology.

Previously, a probe could exhaust hundreds of rays outside a quadrant before
trying any near-reference ray through it. The first outer-sheet torus in the
`tangent-boss` baseline spent over 15M work checking its tube seam. Probes now
search both sides of the quadrant before advancing to finer rays.

`unequal-ordered` reaches checked export at the unchanged 1M request limit.
Raw STEP has one valid solid and ten faces, matches independent OCCT volume
and passes ten material witnesses. The native fillet definition is unchanged.
The default 17 digits and 0.0000001 mm tolerance remain unchanged, as do fitted
fillet bounds, radius, material side and overflow policy.

No public API or recipe field changes. No new fitted carrier, UI work or
SolidDraft pin change. The other consumer-budget cases remain coverage debt;
a faster probe is not general torus or fillet support.
