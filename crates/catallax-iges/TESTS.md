# Coverage inventory

- Public box, bored plate, P2 collar, cylinder and conical frustum: deterministic
  bytes, face counts, shared topology, retained exact supports.
- Native imports: base IGESToBRep without healing; valid solids, exact surface
  kinds, original face counts, two edge uses, shared seams, MM units and signed
  orientation/volume diagnostics. Additional rounded box and rationally posed
  cylinder. Both OCCT 7.8 and 8.0 are exercised independently of the writer.
- Bound controls: nonterminating rational dimensions, independently decoded
  emitted vertices, inclusive acceptance at the rational report bound, rejection
  below it, composition of cone trim and common-scale errors, unchanged sink.
- Refusals: sphere/torus surfaces, uncovered cone angle, zero proof allowance,
  invalid options. No partial bytes on any geometric refusal.
- Retention: STEP bytes and reports unchanged before/after IGES for all three
  requested fixtures; frozen-base CLI output is also compared externally.

Remaining coverage includes general algebraic coordinates, nonterminating
circle-frame components, arbitrary cone half-angles, apex caps, periodic sector
trims and nonaligned ring vertices. These are fail-closed limitations; no test
asserts full support for those families. Native read-back is an explicit separate
command, not a skipped success in the default Rust suite.
