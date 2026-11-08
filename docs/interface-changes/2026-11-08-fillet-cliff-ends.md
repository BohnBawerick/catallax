## 2026-11-08: a wider fillet takes the cliff at a narrower earlier round

Pieces: features, step, app.

Change: under roll-over, a straight strip between planes whose foot reaches an earlier round of its convexity, narrower than its ball and tangent to the foot's face, now builds by the cliff rule of recipe section 4.5.1. The strip keeps its planes and runs on, and the faces it meets cut it. The cut walks from that foot across those faces to the other rail's end, along the strip cylinder's intersections with them. Edges inside the removed region go and edges leaving it end on the cut. The strip face keeps `blend.edge`; every cut face keeps its origins, and no `corner.vertex` face is made. Convex and concave strips both build.

- At a mitre of two such rounds, the strip refused with `fillet.unsupported-corner`. The cut now crosses the mitre, which keeps its outer part. Examples are a block's uprights at radius 2 after its top loop at 1, and a pocket's uprights at 3/2 after its floor loop at 1.
- A tangent chain stops at such a round instead of continuing over it, and its straight strip takes the cliff. It refused with `fillet.radius-too-large` and the round's radius as the limit. Examples are a block's top edge at radius 3 after its uprights at 1, a pocket's floor edge at 3 after its uprights at 1, and a top edge at 2 reaching an earlier edge round of 1.

Why: a ball wider than a round cannot roll onto it, so section 4.5.1 defines that end as a cliff. The result is exact: planes, cylinders and their intersection curves.

Migration: SolidDraft sees bodies where it saw those refusals, with no call changes. An arc of a narrower round chosen by name still refuses with `fillet.radius-too-large`. A ball as wide as the rounds still closes with the sphere sector or sphere strip. A ball narrower than a mitre's rounds, which would roll over the mitre to its apex, still refuses with `fillet.unsupported-corner`. The field compiler is unchanged and keeps refusing roll-over across earlier fillets. STEP writes each cut curve between two curved faces as an approximated curve under the existing file bound.
