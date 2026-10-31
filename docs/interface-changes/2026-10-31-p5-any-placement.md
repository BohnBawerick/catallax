## 2026-10-31: P5 junction continuation in any exact placement

Pieces: measure, features
Change: the P5 rolling-junction continuation now answers when the definition is rigidly placed in an oblique frame, instead of refusing with `fillet.unsupported-overflow`. `catallax_measure::face_heights` gives exact values that enclose a height along one direction over a plane, cylinder or cone face bounded by lines and circles, and the continuation uses them to prove a neighbour lies below the mouth plane or above the ball. World boxes stay as the fallback for other faces.
Why: a world box of a tilted mouth circle reaches corners the circle never reaches, so the guard refused faces that were clear of the strip. Dropping the guard would have accepted a tube too short to carry the rail, which still refuses.
Migration: a caller that matched the rotated P5 refusal now sees the manifest/4 record of the [certified exchange](2026-10-31-p5-certified-exchange.md). Rail clear of or tangent to a mouth, and radii past the contact regime, keep their refusals in every frame. No other recipe changes.
