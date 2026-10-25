## 2026-10-22: share zero-system query scheduling without reordering bounds

Piece: 2
Change: curve membership now uses the same separation-first scheduling as point equality. A fixed-size pending array retains each equation's original index, so reordered decisions still feed the original equation units into distance bounds.
Why: a coincident first curve equation could enter expensive elimination before another equation could establish separation.
Migration: no public signatures or default policy changes. Cached signs and metered enclosures select the query order only; the decider still certifies every conclusion. Existing analytic-curve, intersection-carrier, zero-work bound, undefined-coordinate and deadline tests pass.

