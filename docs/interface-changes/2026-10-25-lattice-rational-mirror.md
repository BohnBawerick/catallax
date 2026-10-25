## 2026-10-25: lattice definitions accept rational reflections

Pieces: field lattice leaf
Change: `Lattice::transform` composes an orientation-reversing rational rigid motion into the frame instead of refusing with `field.lattice-frame`. Non-rational motions keep that refusal.
Why: a mirror of a constant or linear lattice is the reflected point set; the frame matrix stays orthogonal, so its transpose is still the inverse.
Migration: none. Mirrored lattice bodies that previously refused now evaluate; callers that matched on the refusal for reflections must drop that case.
