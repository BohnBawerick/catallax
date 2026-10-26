## 2026-10-21: enclose certified root intervals without re-proving their endpoints

Piece: 1
Change: algebraic enclosure checks the certified interval width before evaluating the defining polynomial at an endpoint. That sign is computed only when bisection is required. Refining calls keep the same bisection sequence; already-narrow calls return the same midpoint/radius bounds without an unused polynomial evaluation.
Why: the residual native profiles spent 48/99 and 25/99 samples in algebraic enclosure. Deep isolating intervals can have thousands of denominator bits even when the caller asks for a modest working precision. Re-evaluating their large defining polynomials consumed time and could spend an enclosure budget without improving its bounds.
Migration: no signatures, default budgets or rounding rules change. The regression first isolates `sqrt(2)/(2^607-1)`, then reaches width `1e-220` under a separate 10,000-work policy; independent GMP signed-square inequalities certify the returned bounds.

