## 2026-10-25: retain complete planar conic incidence

Pieces: 2, 5
Change: a full circle or ellipse retained on a plane and a nondegenerate
quadric reuses those support facts when the query point lies on both supports.
An append-only support index follows the recorded definitions through store
clones. Sphere, circular/elliptic cylinder and cone carriers qualify.
Why: the quadric's restriction to the conic plane is a nonzero multiple of the
conic equation, so the retained incidences already prove membership.
Migration: no signatures, codes, allowances or tolerances change. Plane-only
evidence is insufficient. A torus may have two circle sections in one plane
and does not qualify; intersection-curve component checks remain unchanged.
This extends [partial support reuse](2026-10-25-retained-contact-support-facts.md).
