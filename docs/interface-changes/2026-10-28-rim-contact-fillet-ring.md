## 2026-10-28: fillet ring continues over countersink mouths

Pieces: 2 (catallax-shapes), 4 (catallax-features), and consumers of fillet results
Change:
- **Fillet ring.** A closed concave fillet ring between a plane and a perpendicular cylinder now continues over countersink mouths that the rolling ball reaches. Where it previously refused with `fillet.unsupported-overflow`, it now returns:
  - elementary torus pieces between the holes;
  - one `SplineSupport::RollingEnvelope` face per hole;
  - one merged plane loop and split mouth edges.
- **New refusals.** Contacts outside this stratum keep their refusal, or refuse with `shapes.rolling-envelope` when the ball would leave the mouth edge.
- **Shapes.** `RollingEnvelope::reframed` reuses one fit for a patterned hole. The envelope's middle row now lies in the cylinder's tangent plane at the rail, so the carrier degree is (8, 2). The fit also proves exact clearance from the cylinder, the mouth plane and the neighbouring torus pieces.

Why: the P5 flanged housing's junction fillet reaches each countersink mouth. Its original and edited radii now build a valid closed body instead of refusing.
Migration: code that matched on P5 refusing at its fillet now receives a body with four `RollingEnvelope` faces. Display, field, measurement and STEP consumers that do not yet handle that family must refuse it rather than treat it as a torus. Every face of the ring keeps the single blend origin, under the split rule.
