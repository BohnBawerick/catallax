## 2026-10-20: an unaffordable exact attempt does not end a decision

Pieces: 1 (and every piece that decides through it)
Change: piece 1's rule 4 now says that an exact attempt (annihilating polynomial, separation bound, algebraic identity) whose proof would pass `max_work` is skipped while refinement continues. A decision ends `NotExact` only when the next refinement step would pass a limit, so a non-zero value whose sign refinement within both limits can show is always `Proven`. No signature, code or budget default changes.
Why: on main the sign loop stopped with `numbers.budget-exhausted` at 64 to 128 bits when a separation-bound attempt was too expensive. Piece 5's case 57 query at scale `10^-30` (case 13 bite) was refused although its value, about `10^-104` against coordinates near `10^5`, needs about 400 bits. With the attempt skipped instead, the query answers `Proven(Inside)` in 1.5 s.
Migration: piece 1 keeps refining after a failed or unaffordable exact attempt. Callers need no change; some `NotExact` outcomes become `Proven`.

