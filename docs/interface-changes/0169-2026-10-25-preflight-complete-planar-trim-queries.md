## 2026-10-25: preflight complete planar trim queries

Pieces: 5

Change: enable native circular-arc winding only when every trim in the complete
face or signature batch has a native line/circle carrier. Mixed or implicit trim
sets retain the earlier polygon/full-ring proof order and generic fallback.

Why: a partial arc attempt could spend substantial algebraic work before finding
an unsupported later trim, then discard the result. The complete-query carrier
check prevents that speculative work without adding numerical predicates.

Migration: none. Public APIs, work allowances, tolerances, ordered geometry,
names and STEP serialization are unchanged.

