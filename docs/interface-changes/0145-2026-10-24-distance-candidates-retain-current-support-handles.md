## 2026-10-24: distance candidates retain current support handles

Pieces: field and comparison
Change: the distance-feature visitor uses the current-definition `(RPoint3, FaceId)` callback contract, matching `Plan::visit_candidates`. Native draft-region bounds remain available to the query-only sampler. Constructing a distance feature from a draft region still requires a retained trim proof; its constraints alone are not assumed to define a convex intersection.
Why: an origin can survive on several replacement supports, and a retained region can contain multiple material cells.
Migration: certify the emitted handle directly and read its origin from the field registry. No compatibility alias, recipe format, tolerance or work allowance changes.
