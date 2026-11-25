Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-21: clash or clear where the overlap Boolean refuses

Pieces: application queries (`catallax_app::queries`)
Change: where the overlap Boolean refuses, or leaves the overlap enclosure unresolved, `body_relation` decides the pair from boundary samples, grid points and a proof that the two boundaries are apart, following SMLib's `SmApiBrepRelationship` and OpenCASCADE's `BRepExtrema_DistShapeShape::SolidTreatment`. A clash proven from a sample has an overlap enclosure from 4 r^3, for a ball of radius r proven inside both solids, to the volume of the two proven boxes' overlap, and carries the `measure.accuracy-not-reached` mark. A solid wholly inside the other overlaps by its own volume. A clear pair whose faces were proven apart has a gap at least the lower bound that proof found. A pair still undecided keeps the Boolean's refusal code, and its detail ends with `; arrangement:` and what the search could not prove.
Why: roadmap 3.3 step 2. Overlapping copies of P2, P3 and a spline loft, and turned cone, cylinder and torus pairs, refused in the Boolean although a clash or a gap was provable.
Migration: callers that treated every `BodyRelation::Clash` overlap as meeting the requested accuracy must read the answer's marks. Every pair that answered before gives the same answer and work.
