## 2026-10-21: bounded exact root-certificate reuse

Piece: 1
Change: each decider lazily holds at most 128 nonlinear root certificates in a FIFO memo. A hit requires exact equality of every rational polynomial coefficient. Complete roots, multiplicities and isolating intervals are reused unchanged; hash collisions never identify polynomials. Shared polynomial keys avoid copying their coefficients into both lookup and eviction indexes. Reuse charges the original certificate work and retains deadline checks.
Why: a full native cylinder case profile spent 35.75 seconds in 3,286 root calls for 63 polynomial fingerprints; 34.78 seconds repeated successful fingerprints. Re-isolating the same exact polynomial dominated the remaining deadline rather than needing more precision or a larger budget.
Migration: no public signatures, default budgets, repeated-decision work guarantees or geometric tolerances change. The memo is bounded and decider-local. Independent GMP certificates cover a working set larger than the memo, and a warm-cache test retains work identity and expired-deadline refusal.

