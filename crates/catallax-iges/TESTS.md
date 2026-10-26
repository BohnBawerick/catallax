# Coverage inventory

- Public box, bored plate, P2 collar, cylinder and conical frustum: deterministic
  bytes, face counts, shared topology, retained exact supports.
- Public sketch/revolve cylinder and cone 90/180/270-degree sectors: deterministic
  shared topology without extra seams or vertices. Opposite cone axes and inner
  reversed faces; exact recharting through a seam while reversing all edges and
  circle reference rays; independently decoded finite, closed pcurve rectangles.
- Native imports: base IGESToBRep without healing; valid solids, exact surface
  kinds, original face counts, two edge uses, shared seams, MM units and signed
  orientation/volume diagnostics. Additional rounded box and rationally posed
  cylinder. Both OCCT 7.8 and 8.0 are exercised independently of the writer.
  Sector imports also check finite angular/slant domains, matching native curve
  ranges and sampled pcurve orientation/compatibility. Sampling and signed volume
  remain numerical diagnostics, separate from the whole-domain certificate.
- Bound controls: nonterminating rational dimensions, independently decoded
  emitted vertices, inclusive acceptance at the rational report bound, rejection
  below it, composition of cone trim and common-scale errors, unchanged sink.
  Scaled finite cones explicitly add angular, slant and common-scale errors and
  exercise inclusive acceptance and refusal below the reported bound.
- Refusals: sphere/torus surfaces, uncovered cone angle, zero proof allowance,
  invalid options. No partial bytes on any geometric refusal.
  Finite 45-degree algebraic endpoints, uncovered sector cone slopes and apex
  sectors retain `UnrepresentableDefinition` and preserve a prefilled sink.
- Retention: STEP bytes and reports unchanged before/after IGES for all three
  requested fixtures; frozen-base CLI output is also compared externally.

Remaining coverage includes general algebraic coordinates, nonterminating
circle-frame components, arbitrary cone half-angles, apex caps, noncardinal or
nonrectangular periodic sector trims and nonaligned ring vertices. These are fail-closed limitations; no test
asserts full support for those families. Native read-back is an explicit separate
command, not a skipped success in the default Rust suite.
