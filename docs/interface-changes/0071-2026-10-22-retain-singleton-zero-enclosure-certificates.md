## 2026-10-22: retain singleton-zero enclosure certificates

Piece: 1
Change: a known-defined immutable zero singleton, or a zero singleton retained after checked refinement, proves coincidence without algebraic elimination or a nonzero-root separation bound. The existing algebraic proof category is used.
Why: sign decisions discarded exact zero enclosures and could exhaust their allowance reconstructing a zero already certified by interval arithmetic.
Migration: no public signatures, default allowances or deadline changes. A ball merely containing zero proves nothing. Definedness must be established before using a singleton, including when zero multiplies an expression with a hidden undefined divisor. The regression fails before this change at one work unit; empty allowances still stop, and tiny nonzero enclosures remain non-coincident.

