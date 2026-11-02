## 2026-11-02: shell around separated through bores

Pieces: catallax-features, catallax-app, SolidDraft's Catallax adapter.

Change: ordinary B-rep `shell` accepts a convex planar outer envelope with parallel cylindrical through bores between the same two perpendicular caps. Openings may select either or both caps, or leave the shell closed. New inner cylinders use the existing `inner.of` origin form. No recipe version, public type, numerical policy or work allowance changes.

Why: a tray with a drilled through hole previously left the convex-planar shell family even when its cavity has exact planar and cylindrical walls. The construction reads the final input topology and reuses the draft route's complete-ring bore inventory and exact clearance proofs. It does not infer a shell from a Boolean recipe's operation order.

Coverage: expanded bore cylinders must stay strictly separated from each other and from the inset outer walls. A consumed planar cavity keeps the input. Opening bore walls or outer side faces, blind or countersunk holes, fillets, nonconvex outer envelopes and meeting offset envelopes remain unsupported. Adjacency-filtered openings and nearest selection among multiple candidates refuse rather than resolving against filled cap polygons. The independent field compiler is unchanged and still refuses this family.

Migration: existing native `Operation::Shell` and recipe/3 or recipe/4 requests need no schema change. SolidDraft must pin this Catallax revision and continue sending its resolved cap `FaceRef` and positive thickness through the existing shell operation. Keep field comparison optional and display the existing typed refusals outside the documented family. Application UI acceptance, including edit/save/reopen, remains a SolidDraft integration gate rather than a claim of this kernel slice.
