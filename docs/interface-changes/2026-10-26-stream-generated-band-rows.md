## 2026-10-26: stream generated band rows through certification

Pieces: display mesh, application, two-kernel comparison

Change: closed cone and torus bands sample and certify each generated row in one angular merge pass. Each newly sampled endpoint checkpoints a bounded cell instead of checkpointing both triangles again in a separate strip traversal. Extra lower-rail advances, copied seam endpoints and lower tails retain explicit checkpoints. The final boundary strip is unchanged.

Why: the separate intermediate-strip traversal duplicated checkpoints already covered by row sampling. On the public P3 comparison mesh, deterministic work falls from 696,107 to 680,797. Triangles, rounded positions, shared rails, exact face bounds and OBJ bytes remain identical. Every triangle still runs the existing orientation, chord, endpoint-residual, length, quality and angular checks, including exact fallbacks.

Migration: no API, format, requested tolerance, certificate or allowance changes. Consumers may observe lower deterministic work and fewer resource refusals. No field-display or Boolean code changes. This refines the [layered-band implementation](0173-2026-10-26-triangulate-closed-cone-and-torus-bands-in-layers.md) without changing its mesh output.
