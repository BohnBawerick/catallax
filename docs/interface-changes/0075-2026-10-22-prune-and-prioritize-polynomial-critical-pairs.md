## 2026-10-22: prune and prioritize polynomial critical pairs

Piece: 4
Change: Buchberger reduction uses the chain criterion only for pairs with previously established standard representations. Remaining pairs are selected by leading-monomial LCM degree, then monomial and generator indices.
Why: redundant FIFO reductions caused coefficient swell in posed finite curve/surface intersections.
Migration: coefficient operations and domains are unchanged; the criterion uses only leading monomials and resolved-pair identities. Reduced monic bases, root ordering and stored geometry definitions retain the same contracts. A cyclic three-equation regression checks both exact real solutions under all six generator orders; the complete scale and placement family passes.

