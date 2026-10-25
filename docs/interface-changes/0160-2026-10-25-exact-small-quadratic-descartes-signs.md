## 2026-10-25: exact small-quadratic Descartes signs

Pieces: 1

Change: complete rational-search branches can evaluate the three transformed
quadratic coefficients with exact native integers instead of allocating
generic interval Taylor shifts. Admission bounds every input to 20 bits and
every intermediate to signed 63 bits. Other coefficients, fractional endpoints,
zero working-precision caps and policy-bounded reconstruction retain the
existing filtered/exact isolator.

Why: Descartes needs coefficient signs, not a generic polynomial object, when
the exact three-coefficient transform already fits one-limb arithmetic.

Migration: none. Public signatures, root-selection proofs, tolerance rules and
work allowances are unchanged. Arithmetic work may decrease.

