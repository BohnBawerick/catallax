## 2026-10-24: oriented planar draft components

Pieces: features, field, compare

Change: `draft::arrangement` adds exact skin patches, planar subdivisions and connected material components. B-rep source wires and field sweep definitions independently supply the input skin. Exact patch winding labels the new plane arrangement; retained components determine the final topology. Connected coplanar patches merge their origins, while disconnected coplanar patches retain distinct source names. Removed-face history is computed across all retained components, so a surviving split face is not reported as consumed.

The field preserves a true-distance certificate for unions of the separated components established by draft. Ordinary Boolean unions remain residuals. Compaction and rigid copies retain the certificate. Native samples and distance queries use surviving trimmed patches, not internal subdivision seams.

Coverage added: nonconvex line-bounded planar skins whose incident replacement carriers have consistent finite intersections. Non-simple or singular moved vertex fans still report explicit unfinished coverage; convex planar cells retain their direct support-intersection path. Mixed curved and Boolean-input field arrangements remain unfinished. This entry corrects the narrower planar coverage recorded above; it does not claim the full P1 latency exit.

Prior art: [OCCT draft modification](https://occt3d.com/dev/doc/refman/html/class_draft___modification.html) separates replacement surfaces, curves and points, but [its public topology-preserving restriction](https://occt3d.com/dev/doc/refman/html/class_b_rep_offset_a_p_i___draft_angle.html) cannot implement swallowed faces. [Shapiro and Vossler's arrangement treatment](https://spatial.engr.wisc.edu/wp-content/uploads/sites/715/2014/04/0000-4.pdf) motivates explicit cell labels. Patch-local provenance avoids assigning every coplanar source name to every output patch. No code is copied.

Migration: no recipe version, public callback signature, tolerance or budget changes. Consumers of `Value` can receive a certified `TrueDistance` for a separated draft-component union that previously had only residual semantics.

