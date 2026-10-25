## 2026-10-22: scale preliminary coefficient enclosures

Piece: 5
Change: the probe-coefficient sign filter first tries the cached conservative interval width when it is looser than the refinement target. Further refinement is bounded relative to that width as well as by the existing absolute target.
Why: at large geometric scales, an absolute width alone forced unnecessary precision before testing a coefficient's sign.
Migration: this changes only a preliminary metered enclosure request. A zero-containing ball still proves nothing, and the exact fallback is unchanged. No geometric tolerance, default allowance, public signature or deadline changes.

