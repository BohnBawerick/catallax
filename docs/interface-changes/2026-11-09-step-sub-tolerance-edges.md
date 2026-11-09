## 2026-11-09: collapse or refuse short straight STEP edges

Piece: 7

STEP export now refuses with `step.sub-tolerance-edge` when a written straight
edge at or below the 1e-7 mm reader floor remains after export repairs. It writes
no bytes on refusal. Construction, exact model topology, tolerance policy and
measurement stay unchanged.

A redundant boundary subdivision can instead collapse without changing the
boundary. The subdivision vertex must have exactly two straight incident edges
on the same carrier, with endpoints in order along the line. Both affected wires
must remain ordinary closed wires, and the copied body must pass `Model::check`.
The regularization record reports `removed_faces: 0`, `removed_edges: 1` and zero
deviation. Its `merged_faces` name the two affected faces without joining them.

The exact witnesses are in `crates/catallax-step/examples/sub_tolerance.rs`:

- A 10 by 10 by 10 mm block with its first boundary edge split at the fraction
  `3.5e-11` of its length. The short segment is 3.5e-10 mm. OCCT 8.0.1.1.0 rejects
  the original raw translation; its standard reader repairs it. Export now
  removes that redundant vertex and writes the original cube boundary.
- A `3.5e-10` by 10 by 10 mm block. Both raw and standard OCCT readers reject the
  original export. Collapsing its short edges would remove its volume, so export
  refuses it.

The original files, generated before the fix, are retained under
`tools/step-regressions/fixtures/sub-tolerance/`. The OCCT regression checks their
raw invalidity before checking the current exporter. The two tests failed before
the implementation changed.

Migration: a uniformly tiny body that previously exported unreadable straight
edges now refuses. Callers must handle the named refusal, as they handle other
export refusals. A looser file tolerance does not lower the reader floor.
