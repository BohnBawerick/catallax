# Algebraic revolutions

The B-rep revolution constructor accepts rational-degree angles in `(0, 360]`
through the shared exact sine/cosine implementation. Rotated vertices, cap
planes, profile lines and profile circles retain algebraic coordinates. Axis
normalization and generated spherical and ring-toroidal radii also remain exact.
The original analytic surface kinds are retained; no fitted or sampled surface
replaces the definition.

Rational radii and rational full/cardinal frames retain their existing
representation. Ring-torus coverage is checked by comparing the rational
squared major and minor radii. Horn and spindle sheets still return
`features.unsupported-spindle-torus`. Invalid angles and axes keep their existing
recipe errors.

Public examples are in
`crates/catallax-features/tests/fixtures/algebraic-revolutions`. The focused
`revolutions` and `algebraic_revolutions` integration tests cover analytic mass
controls, exact angular membership, topology, origins, parameter-edit names,
default requested mass properties, certified B-rep display, and STEP export.
Measurement and display remain explicit consumer requests; ordinary construction
does not collect qualification measurements or timing information.

## Pending field integration

This draft includes the independent field angle implementation in
`crates/catallax-field/src/features/revolution.rs`. Its shared dispatch hook is
pending integration. Until that hook lands, the field compiler retains its
previous algebraic revolution refusal codes.

The native field display charts have a further coverage restriction: they store
rational chart points and witnesses. A general algebraic cap plane cannot use
that representation. Enabling exact field construction and angular queries does
not by itself establish native field display coverage. That consumer needs an
algebraic chart/witness extension before this slice can be declared complete.

Rational-only comparison candidate APIs also need an explicit coverage refusal
for algebraic cap samples; rounding a rotated point is not an exact sample.
Independent STEP import establishes exchange validity, not a certified volume
or deviation enclosure. Public regression replay does not replace qualification
against an unavailable external release corpus.
